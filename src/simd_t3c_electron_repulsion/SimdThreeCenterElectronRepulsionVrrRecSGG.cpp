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


#include "SimdThreeCenterElectronRepulsionVrrRecSGG.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_sgg_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sfg0,
                                                          const size_t sff, const size_t sfg1,
                                                          const size_t sgd0, const size_t sgd1,
                                                          const size_t sgf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 1.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfg0_0 = buffer.data(sfg0 + 0);
    const auto *sfg0_3 = buffer.data(sfg0 + 3);
    const auto *sfg0_5 = buffer.data(sfg0 + 5);
    const auto *sfg0_10 = buffer.data(sfg0 + 10);
    const auto *sfg0_14 = buffer.data(sfg0 + 14);
    const auto *sfg0_18 = buffer.data(sfg0 + 18);
    const auto *sfg0_25 = buffer.data(sfg0 + 25);
    const auto *sfg0_30 = buffer.data(sfg0 + 30);
    const auto *sfg0_35 = buffer.data(sfg0 + 35);
    const auto *sfg0_44 = buffer.data(sfg0 + 44);
    const auto *sfg0_45 = buffer.data(sfg0 + 45);
    const auto *sfg0_48 = buffer.data(sfg0 + 48);
    const auto *sfg0_75 = buffer.data(sfg0 + 75);
    const auto *sfg0_90 = buffer.data(sfg0 + 90);
    const auto *sfg0_93 = buffer.data(sfg0 + 93);
    const auto *sfg0_95 = buffer.data(sfg0 + 95);
    const auto *sfg0_100 = buffer.data(sfg0 + 100);
    const auto *sfg0_102 = buffer.data(sfg0 + 102);
    const auto *sfg0_104 = buffer.data(sfg0 + 104);
    const auto *sfg0_110 = buffer.data(sfg0 + 110);
    const auto *sfg0_115 = buffer.data(sfg0 + 115);
    const auto *sfg0_117 = buffer.data(sfg0 + 117);
    const auto *sfg0_119 = buffer.data(sfg0 + 119);
    const auto *sfg0_123 = buffer.data(sfg0 + 123);

    const auto *sff_0 = buffer.data(sff + 0);
    const auto *sff_1 = buffer.data(sff + 1);
    const auto *sff_2 = buffer.data(sff + 2);
    const auto *sff_3 = buffer.data(sff + 3);
    const auto *sff_5 = buffer.data(sff + 5);
    const auto *sff_6 = buffer.data(sff + 6);
    const auto *sff_7 = buffer.data(sff + 7);
    const auto *sff_8 = buffer.data(sff + 8);
    const auto *sff_9 = buffer.data(sff + 9);
    const auto *sff_10 = buffer.data(sff + 10);
    const auto *sff_12 = buffer.data(sff + 12);
    const auto *sff_16 = buffer.data(sff + 16);
    const auto *sff_17 = buffer.data(sff + 17);
    const auto *sff_18 = buffer.data(sff + 18);
    const auto *sff_19 = buffer.data(sff + 19);
    const auto *sff_20 = buffer.data(sff + 20);
    const auto *sff_22 = buffer.data(sff + 22);
    const auto *sff_26 = buffer.data(sff + 26);
    const auto *sff_27 = buffer.data(sff + 27);
    const auto *sff_28 = buffer.data(sff + 28);
    const auto *sff_29 = buffer.data(sff + 29);
    const auto *sff_30 = buffer.data(sff + 30);
    const auto *sff_32 = buffer.data(sff + 32);
    const auto *sff_33 = buffer.data(sff + 33);
    const auto *sff_35 = buffer.data(sff + 35);
    const auto *sff_36 = buffer.data(sff + 36);
    const auto *sff_37 = buffer.data(sff + 37);
    const auto *sff_38 = buffer.data(sff + 38);
    const auto *sff_39 = buffer.data(sff + 39);
    const auto *sff_40 = buffer.data(sff + 40);
    const auto *sff_42 = buffer.data(sff + 42);
    const auto *sff_46 = buffer.data(sff + 46);
    const auto *sff_47 = buffer.data(sff + 47);
    const auto *sff_48 = buffer.data(sff + 48);
    const auto *sff_49 = buffer.data(sff + 49);
    const auto *sff_50 = buffer.data(sff + 50);
    const auto *sff_52 = buffer.data(sff + 52);
    const auto *sff_53 = buffer.data(sff + 53);
    const auto *sff_55 = buffer.data(sff + 55);
    const auto *sff_56 = buffer.data(sff + 56);
    const auto *sff_57 = buffer.data(sff + 57);
    const auto *sff_58 = buffer.data(sff + 58);
    const auto *sff_59 = buffer.data(sff + 59);
    const auto *sff_60 = buffer.data(sff + 60);
    const auto *sff_63 = buffer.data(sff + 63);
    const auto *sff_65 = buffer.data(sff + 65);
    const auto *sff_66 = buffer.data(sff + 66);
    const auto *sff_67 = buffer.data(sff + 67);
    const auto *sff_68 = buffer.data(sff + 68);
    const auto *sff_69 = buffer.data(sff + 69);
    const auto *sff_75 = buffer.data(sff + 75);
    const auto *sff_76 = buffer.data(sff + 76);
    const auto *sff_77 = buffer.data(sff + 77);
    const auto *sff_78 = buffer.data(sff + 78);
    const auto *sff_79 = buffer.data(sff + 79);
    const auto *sff_83 = buffer.data(sff + 83);

    const auto *sfg1_0 = buffer.data(sfg1 + 0);
    const auto *sfg1_3 = buffer.data(sfg1 + 3);
    const auto *sfg1_5 = buffer.data(sfg1 + 5);
    const auto *sfg1_10 = buffer.data(sfg1 + 10);
    const auto *sfg1_14 = buffer.data(sfg1 + 14);
    const auto *sfg1_18 = buffer.data(sfg1 + 18);
    const auto *sfg1_25 = buffer.data(sfg1 + 25);
    const auto *sfg1_30 = buffer.data(sfg1 + 30);
    const auto *sfg1_35 = buffer.data(sfg1 + 35);
    const auto *sfg1_44 = buffer.data(sfg1 + 44);
    const auto *sfg1_45 = buffer.data(sfg1 + 45);
    const auto *sfg1_48 = buffer.data(sfg1 + 48);
    const auto *sfg1_75 = buffer.data(sfg1 + 75);
    const auto *sfg1_90 = buffer.data(sfg1 + 90);
    const auto *sfg1_93 = buffer.data(sfg1 + 93);
    const auto *sfg1_95 = buffer.data(sfg1 + 95);
    const auto *sfg1_100 = buffer.data(sfg1 + 100);
    const auto *sfg1_102 = buffer.data(sfg1 + 102);
    const auto *sfg1_104 = buffer.data(sfg1 + 104);
    const auto *sfg1_110 = buffer.data(sfg1 + 110);
    const auto *sfg1_115 = buffer.data(sfg1 + 115);
    const auto *sfg1_117 = buffer.data(sfg1 + 117);
    const auto *sfg1_119 = buffer.data(sfg1 + 119);
    const auto *sfg1_123 = buffer.data(sfg1 + 123);

    const auto *sgd0_0 = buffer.data(sgd0 + 0);
    const auto *sgd0_3 = buffer.data(sgd0 + 3);
    const auto *sgd0_5 = buffer.data(sgd0 + 5);
    const auto *sgd0_9 = buffer.data(sgd0 + 9);
    const auto *sgd0_11 = buffer.data(sgd0 + 11);
    const auto *sgd0_17 = buffer.data(sgd0 + 17);
    const auto *sgd0_18 = buffer.data(sgd0 + 18);
    const auto *sgd0_21 = buffer.data(sgd0 + 21);
    const auto *sgd0_23 = buffer.data(sgd0 + 23);
    const auto *sgd0_29 = buffer.data(sgd0 + 29);
    const auto *sgd0_30 = buffer.data(sgd0 + 30);
    const auto *sgd0_33 = buffer.data(sgd0 + 33);
    const auto *sgd0_35 = buffer.data(sgd0 + 35);

    const auto *sgd1_0 = buffer.data(sgd1 + 0);
    const auto *sgd1_3 = buffer.data(sgd1 + 3);
    const auto *sgd1_5 = buffer.data(sgd1 + 5);
    const auto *sgd1_9 = buffer.data(sgd1 + 9);
    const auto *sgd1_11 = buffer.data(sgd1 + 11);
    const auto *sgd1_17 = buffer.data(sgd1 + 17);
    const auto *sgd1_18 = buffer.data(sgd1 + 18);
    const auto *sgd1_21 = buffer.data(sgd1 + 21);
    const auto *sgd1_23 = buffer.data(sgd1 + 23);
    const auto *sgd1_29 = buffer.data(sgd1 + 29);
    const auto *sgd1_30 = buffer.data(sgd1 + 30);
    const auto *sgd1_33 = buffer.data(sgd1 + 33);
    const auto *sgd1_35 = buffer.data(sgd1 + 35);

    const auto *sgf_0 = buffer.data(sgf + 0);
    const auto *sgf_2 = buffer.data(sgf + 2);
    const auto *sgf_3 = buffer.data(sgf + 3);
    const auto *sgf_5 = buffer.data(sgf + 5);
    const auto *sgf_6 = buffer.data(sgf + 6);
    const auto *sgf_7 = buffer.data(sgf + 7);
    const auto *sgf_8 = buffer.data(sgf + 8);
    const auto *sgf_9 = buffer.data(sgf + 9);
    const auto *sgf_10 = buffer.data(sgf + 10);
    const auto *sgf_12 = buffer.data(sgf + 12);
    const auto *sgf_16 = buffer.data(sgf + 16);
    const auto *sgf_17 = buffer.data(sgf + 17);
    const auto *sgf_18 = buffer.data(sgf + 18);
    const auto *sgf_19 = buffer.data(sgf + 19);
    const auto *sgf_20 = buffer.data(sgf + 20);
    const auto *sgf_22 = buffer.data(sgf + 22);
    const auto *sgf_26 = buffer.data(sgf + 26);
    const auto *sgf_27 = buffer.data(sgf + 27);
    const auto *sgf_28 = buffer.data(sgf + 28);
    const auto *sgf_29 = buffer.data(sgf + 29);
    const auto *sgf_30 = buffer.data(sgf + 30);
    const auto *sgf_32 = buffer.data(sgf + 32);
    const auto *sgf_33 = buffer.data(sgf + 33);
    const auto *sgf_35 = buffer.data(sgf + 35);
    const auto *sgf_36 = buffer.data(sgf + 36);
    const auto *sgf_37 = buffer.data(sgf + 37);
    const auto *sgf_38 = buffer.data(sgf + 38);
    const auto *sgf_39 = buffer.data(sgf + 39);
    const auto *sgf_40 = buffer.data(sgf + 40);
    const auto *sgf_42 = buffer.data(sgf + 42);
    const auto *sgf_46 = buffer.data(sgf + 46);
    const auto *sgf_47 = buffer.data(sgf + 47);
    const auto *sgf_48 = buffer.data(sgf + 48);
    const auto *sgf_49 = buffer.data(sgf + 49);
    const auto *sgf_50 = buffer.data(sgf + 50);
    const auto *sgf_52 = buffer.data(sgf + 52);
    const auto *sgf_53 = buffer.data(sgf + 53);
    const auto *sgf_55 = buffer.data(sgf + 55);
    const auto *sgf_56 = buffer.data(sgf + 56);
    const auto *sgf_57 = buffer.data(sgf + 57);
    const auto *sgf_58 = buffer.data(sgf + 58);
    const auto *sgf_59 = buffer.data(sgf + 59);
    const auto *sgf_60 = buffer.data(sgf + 60);
    const auto *sgf_62 = buffer.data(sgf + 62);
    const auto *sgf_66 = buffer.data(sgf + 66);
    const auto *sgf_67 = buffer.data(sgf + 67);
    const auto *sgf_68 = buffer.data(sgf + 68);
    const auto *sgf_69 = buffer.data(sgf + 69);
    const auto *sgf_70 = buffer.data(sgf + 70);
    const auto *sgf_72 = buffer.data(sgf + 72);
    const auto *sgf_76 = buffer.data(sgf + 76);
    const auto *sgf_77 = buffer.data(sgf + 77);
    const auto *sgf_78 = buffer.data(sgf + 78);
    const auto *sgf_79 = buffer.data(sgf + 79);
    const auto *sgf_80 = buffer.data(sgf + 80);
    const auto *sgf_82 = buffer.data(sgf + 82);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, sff_0, sff_3, sgd0_0, sgd0_3, \
                         sgd1_0, sgd1_3, sgf_0, sgf_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sff_0[k]
                 + f_1 * sgd0_0[k]
                 - f_2 * sgd1_0[k]
                 + f_3 * pc_x[k] * sgf_0[k];

        t_1[k] = f_3 * pc_y[k] * sgf_0[k];

        t_2[k] = f_3 * pc_z[k] * sgf_0[k];

        t_3[k] = f_0 * sff_3[k]
                 + f_4 * sgd0_3[k]
                 - f_5 * sgd1_3[k]
                 + f_3 * pc_x[k] * sgf_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pc_x, pc_y, sff_5, sff_6, sff_7, sgd0_5, sgd1_5, \
                         sgf_2, sgf_5, sgf_6, sgf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * sgf_2[k];

        t_5[k] = f_0 * sff_5[k]
                 + f_4 * sgd0_5[k]
                 - f_5 * sgd1_5[k]
                 + f_3 * pc_x[k] * sgf_5[k];

        t_6[k] = f_0 * sff_6[k]
                 + f_3 * pc_x[k] * sgf_6[k];

        t_7[k] = f_0 * sff_7[k]
                 + f_3 * pc_x[k] * sgf_7[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pc_x, pc_y, pc_z, sff_8, sff_9, sgd0_3, sgd1_3, \
                         sgf_6, sgf_8, sgf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * sff_8[k]
                 + f_3 * pc_x[k] * sgf_8[k];

        t_9[k] = f_0 * sff_9[k]
                 + f_3 * pc_x[k] * sgf_9[k];

        t_10[k] = f_1 * sgd0_3[k]
                  - f_2 * sgd1_3[k]
                  + f_3 * pc_y[k] * sgf_6[k];

        t_11[k] = f_3 * pc_z[k] * sgf_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, pb_y, pc_y, pc_z, sfg0_0, sff_0, \
                         sfg1_0, sgd0_5, sgd1_5, sgf_8, sgf_9, sgf_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_4 * sgd0_5[k]
                  - f_5 * sgd1_5[k]
                  + f_3 * pc_y[k] * sgf_8[k];

        t_13[k] = f_3 * pc_y[k] * sgf_9[k];

        t_14[k] = f_1 * sgd0_5[k]
                  - f_2 * sgd1_5[k]
                  + f_3 * pc_z[k] * sgf_9[k];

        t_15[k] = pb_y[k] * sfg0_0[k]
                  - f_6 * pc_y[k] * sfg1_0[k];

        t_16[k] = f_7 * sff_0[k]
                  + f_3 * pc_y[k] * sgf_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pb_y, pc_y, pc_z, sfg0_3, sfg0_5, sff_1, \
                         sff_2, sfg1_3, sfg1_5, sgf_10, sgf_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_3 * pc_z[k] * sgf_10[k];

        t_18[k] = pb_y[k] * sfg0_3[k]
                  + f_8 * sff_1[k]
                  - f_6 * pc_y[k] * sfg1_3[k];

        t_19[k] = f_7 * sff_2[k]
                  + f_3 * pc_y[k] * sgf_12[k];

        t_20[k] = pb_y[k] * sfg0_5[k]
                  - f_6 * pc_y[k] * sfg1_5[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_x, sff_16, sff_17, sff_18, sff_19, sgf_16, \
                         sgf_17, sgf_18, sgf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_9 * sff_16[k]
                  + f_3 * pc_x[k] * sgf_16[k];

        t_22[k] = f_9 * sff_17[k]
                  + f_3 * pc_x[k] * sgf_17[k];

        t_23[k] = f_9 * sff_18[k]
                  + f_3 * pc_x[k] * sgf_18[k];

        t_24[k] = f_9 * sff_19[k]
                  + f_3 * pc_x[k] * sgf_19[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pc_y, pc_z, sff_6, sff_8, sff_9, sgd0_9, \
                         sgd0_11, sgd1_9, sgd1_11, sgf_16, sgf_18, \
                         sgf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_7 * sff_6[k]
                  + f_1 * sgd0_9[k]
                  - f_2 * sgd1_9[k]
                  + f_3 * pc_y[k] * sgf_16[k];

        t_26[k] = f_3 * pc_z[k] * sgf_16[k];

        t_27[k] = f_7 * sff_8[k]
                  + f_4 * sgd0_11[k]
                  - f_5 * sgd1_11[k]
                  + f_3 * pc_y[k] * sgf_18[k];

        t_28[k] = f_7 * sff_9[k]
                  + f_3 * pc_y[k] * sgf_19[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pb_y, pb_z, pc_y, pc_z, sfg0_0, sfg0_14, \
                         sff_0, sfg1_0, sfg1_14, sgf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pb_y[k] * sfg0_14[k]
                  - f_6 * pc_y[k] * sfg1_14[k];

        t_30[k] = pb_z[k] * sfg0_0[k]
                  - f_6 * pc_z[k] * sfg1_0[k];

        t_31[k] = f_3 * pc_y[k] * sgf_20[k];

        t_32[k] = f_7 * sff_0[k]
                  + f_3 * pc_z[k] * sgf_20[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pb_z, pc_x, pc_y, pc_z, sfg0_3, sfg0_5, \
                         sff_2, sff_26, sfg1_3, sfg1_5, sgf_22, \
                         sgf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pb_z[k] * sfg0_3[k]
                  - f_6 * pc_z[k] * sfg1_3[k];

        t_34[k] = f_3 * pc_y[k] * sgf_22[k];

        t_35[k] = pb_z[k] * sfg0_5[k]
                  + f_8 * sff_2[k]
                  - f_6 * pc_z[k] * sfg1_5[k];

        t_36[k] = f_9 * sff_26[k]
                  + f_3 * pc_x[k] * sgf_26[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pb_z, pc_x, pc_z, sfg0_10, sff_27, sff_28, \
                         sff_29, sfg1_10, sgf_27, sgf_28, sgf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_9 * sff_27[k]
                  + f_3 * pc_x[k] * sgf_27[k];

        t_38[k] = f_9 * sff_28[k]
                  + f_3 * pc_x[k] * sgf_28[k];

        t_39[k] = f_9 * sff_29[k]
                  + f_3 * pc_x[k] * sgf_29[k];

        t_40[k] = pb_z[k] * sfg0_10[k]
                  - f_6 * pc_z[k] * sfg1_10[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pc_y, pc_z, sff_6, sff_9, sgd0_17, sgd1_17, \
                         sgf_26, sgf_28, sgf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_7 * sff_6[k]
                  + f_3 * pc_z[k] * sgf_26[k];

        t_42[k] = f_4 * sgd0_17[k]
                  - f_5 * sgd1_17[k]
                  + f_3 * pc_y[k] * sgf_28[k];

        t_43[k] = f_3 * pc_y[k] * sgf_29[k];

        t_44[k] = f_7 * sff_9[k]
                  + f_1 * sgd0_17[k]
                  - f_2 * sgd1_17[k]
                  + f_3 * pc_z[k] * sgf_29[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pc_x, pc_y, pc_z, sff_10, sff_30, sff_33, \
                         sgd0_18, sgd0_21, sgd1_18, sgd1_21, sgf_30, \
                         sgf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_8 * sff_30[k]
                  + f_1 * sgd0_18[k]
                  - f_2 * sgd1_18[k]
                  + f_3 * pc_x[k] * sgf_30[k];

        t_46[k] = f_8 * sff_10[k]
                  + f_3 * pc_y[k] * sgf_30[k];

        t_47[k] = f_3 * pc_z[k] * sgf_30[k];

        t_48[k] = f_8 * sff_33[k]
                  + f_4 * sgd0_21[k]
                  - f_5 * sgd1_21[k]
                  + f_3 * pc_x[k] * sgf_33[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pc_x, pc_y, sff_12, sff_35, sff_36, sff_37, \
                         sgd0_23, sgd1_23, sgf_32, sgf_35, sgf_36, \
                         sgf_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_8 * sff_12[k]
                  + f_3 * pc_y[k] * sgf_32[k];

        t_50[k] = f_8 * sff_35[k]
                  + f_4 * sgd0_23[k]
                  - f_5 * sgd1_23[k]
                  + f_3 * pc_x[k] * sgf_35[k];

        t_51[k] = f_8 * sff_36[k]
                  + f_3 * pc_x[k] * sgf_36[k];

        t_52[k] = f_8 * sff_37[k]
                  + f_3 * pc_x[k] * sgf_37[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pc_x, pc_y, pc_z, sff_16, sff_38, sff_39, \
                         sgd0_21, sgd1_21, sgf_36, sgf_38, sgf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_8 * sff_38[k]
                  + f_3 * pc_x[k] * sgf_38[k];

        t_54[k] = f_8 * sff_39[k]
                  + f_3 * pc_x[k] * sgf_39[k];

        t_55[k] = f_8 * sff_16[k]
                  + f_1 * sgd0_21[k]
                  - f_2 * sgd1_21[k]
                  + f_3 * pc_y[k] * sgf_36[k];

        t_56[k] = f_3 * pc_z[k] * sgf_36[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pb_y, pc_y, pc_z, sfg0_30, sff_18, sff_19, \
                         sfg1_30, sgd0_23, sgd1_23, sgf_38, sgf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_8 * sff_18[k]
                  + f_4 * sgd0_23[k]
                  - f_5 * sgd1_23[k]
                  + f_3 * pc_y[k] * sgf_38[k];

        t_58[k] = f_8 * sff_19[k]
                  + f_3 * pc_y[k] * sgf_39[k];

        t_59[k] = f_1 * sgd0_23[k]
                  - f_2 * sgd1_23[k]
                  + f_3 * pc_z[k] * sgf_39[k];

        t_60[k] = pb_y[k] * sfg0_30[k]
                  - f_6 * pc_y[k] * sfg1_30[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pb_z, pc_y, pc_z, sfg0_18, sff_10, sff_20, \
                         sff_22, sfg1_18, sgf_40, sgf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_7 * sff_20[k]
                  + f_3 * pc_y[k] * sgf_40[k];

        t_62[k] = f_7 * sff_10[k]
                  + f_3 * pc_z[k] * sgf_40[k];

        t_63[k] = pb_z[k] * sfg0_18[k]
                  - f_6 * pc_z[k] * sfg1_18[k];

        t_64[k] = f_7 * sff_22[k]
                  + f_3 * pc_y[k] * sgf_42[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pb_y, pc_x, pc_y, sfg0_35, sff_46, sff_47, \
                         sff_48, sfg1_35, sgf_46, sgf_47, sgf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pb_y[k] * sfg0_35[k]
                  - f_6 * pc_y[k] * sfg1_35[k];

        t_66[k] = f_8 * sff_46[k]
                  + f_3 * pc_x[k] * sgf_46[k];

        t_67[k] = f_8 * sff_47[k]
                  + f_3 * pc_x[k] * sgf_47[k];

        t_68[k] = f_8 * sff_48[k]
                  + f_3 * pc_x[k] * sgf_48[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pb_z, pc_x, pc_z, sfg0_25, sff_16, sff_49, sfg1_25, \
                         sgf_46, sgf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_8 * sff_49[k]
                  + f_3 * pc_x[k] * sgf_49[k];

        t_70[k] = pb_z[k] * sfg0_25[k]
                  - f_6 * pc_z[k] * sfg1_25[k];

        t_71[k] = f_7 * sff_16[k]
                  + f_3 * pc_z[k] * sgf_46[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pb_y, pc_y, sfg0_44, sff_28, sff_29, sfg1_44, \
                         sgd0_29, sgd1_29, sgf_48, sgf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_7 * sff_28[k]
                  + f_4 * sgd0_29[k]
                  - f_5 * sgd1_29[k]
                  + f_3 * pc_y[k] * sgf_48[k];

        t_73[k] = f_7 * sff_29[k]
                  + f_3 * pc_y[k] * sgf_49[k];

        t_74[k] = pb_y[k] * sfg0_44[k]
                  - f_6 * pc_y[k] * sfg1_44[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pc_x, pc_y, pc_z, sff_20, sff_50, sff_53, \
                         sgd0_30, sgd0_33, sgd1_30, sgd1_33, sgf_50, \
                         sgf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_8 * sff_50[k]
                  + f_1 * sgd0_30[k]
                  - f_2 * sgd1_30[k]
                  + f_3 * pc_x[k] * sgf_50[k];

        t_76[k] = f_3 * pc_y[k] * sgf_50[k];

        t_77[k] = f_8 * sff_20[k]
                  + f_3 * pc_z[k] * sgf_50[k];

        t_78[k] = f_8 * sff_53[k]
                  + f_4 * sgd0_33[k]
                  - f_5 * sgd1_33[k]
                  + f_3 * pc_x[k] * sgf_53[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pc_x, pc_y, sff_55, sff_56, sff_57, sgd0_35, \
                         sgd1_35, sgf_52, sgf_55, sgf_56, sgf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_3 * pc_y[k] * sgf_52[k];

        t_80[k] = f_8 * sff_55[k]
                  + f_4 * sgd0_35[k]
                  - f_5 * sgd1_35[k]
                  + f_3 * pc_x[k] * sgf_55[k];

        t_81[k] = f_8 * sff_56[k]
                  + f_3 * pc_x[k] * sgf_56[k];

        t_82[k] = f_8 * sff_57[k]
                  + f_3 * pc_x[k] * sgf_57[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pc_x, pc_y, pc_z, sff_26, sff_58, sff_59, \
                         sgd0_33, sgd1_33, sgf_56, sgf_58, sgf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_8 * sff_58[k]
                  + f_3 * pc_x[k] * sgf_58[k];

        t_84[k] = f_8 * sff_59[k]
                  + f_3 * pc_x[k] * sgf_59[k];

        t_85[k] = f_1 * sgd0_33[k]
                  - f_2 * sgd1_33[k]
                  + f_3 * pc_y[k] * sgf_56[k];

        t_86[k] = f_8 * sff_26[k]
                  + f_3 * pc_z[k] * sgf_56[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pb_x, pc_x, pc_y, pc_z, sfg0_90, sff_29, \
                         sff_60, sfg1_90, sgd0_35, sgd1_35, sgf_58, \
                         sgf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_4 * sgd0_35[k]
                  - f_5 * sgd1_35[k]
                  + f_3 * pc_y[k] * sgf_58[k];

        t_88[k] = f_3 * pc_y[k] * sgf_59[k];

        t_89[k] = f_8 * sff_29[k]
                  + f_1 * sgd0_35[k]
                  - f_2 * sgd1_35[k]
                  + f_3 * pc_z[k] * sgf_59[k];

        t_90[k] = pb_x[k] * sfg0_90[k]
                  + f_0 * sff_60[k]
                  - f_6 * pc_x[k] * sfg1_90[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pb_x, pc_x, pc_y, pc_z, sfg0_93, sff_30, \
                         sff_32, sff_63, sfg1_93, sgf_60, sgf_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_9 * sff_30[k]
                  + f_3 * pc_y[k] * sgf_60[k];

        t_92[k] = f_3 * pc_z[k] * sgf_60[k];

        t_93[k] = pb_x[k] * sfg0_93[k]
                  + f_8 * sff_63[k]
                  - f_6 * pc_x[k] * sfg1_93[k];

        t_94[k] = f_9 * sff_32[k]
                  + f_3 * pc_y[k] * sgf_62[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pb_x, pc_x, sfg0_95, sff_65, sff_66, sff_67, \
                         sff_68, sfg1_95, sgf_66, sgf_67, sgf_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = pb_x[k] * sfg0_95[k]
                  + f_8 * sff_65[k]
                  - f_6 * pc_x[k] * sfg1_95[k];

        t_96[k] = f_7 * sff_66[k]
                  + f_3 * pc_x[k] * sgf_66[k];

        t_97[k] = f_7 * sff_67[k]
                  + f_3 * pc_x[k] * sgf_67[k];

        t_98[k] = f_7 * sff_68[k]
                  + f_3 * pc_x[k] * sgf_68[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pb_x, pc_x, pc_z, sfg0_100, sfg0_102, \
                         sff_69, sfg1_100, sfg1_102, sgf_66, sgf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_7 * sff_69[k]
                  + f_3 * pc_x[k] * sgf_69[k];

        t_100[k] = pb_x[k] * sfg0_100[k]
                   - f_6 * pc_x[k] * sfg1_100[k];

        t_101[k] = f_3 * pc_z[k] * sgf_66[k];

        t_102[k] = pb_x[k] * sfg0_102[k]
                   - f_6 * pc_x[k] * sfg1_102[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, pb_x, pb_z, pc_x, pc_y, pc_z, sfg0_45, sfg0_104, \
                         sff_39, sfg1_45, sfg1_104, sgf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_9 * sff_39[k]
                   + f_3 * pc_y[k] * sgf_69[k];

        t_104[k] = pb_x[k] * sfg0_104[k]
                   - f_6 * pc_x[k] * sfg1_104[k];

        t_105[k] = pb_z[k] * sfg0_45[k]
                   - f_6 * pc_z[k] * sfg1_45[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pb_z, pc_y, pc_z, sfg0_48, sff_30, \
                         sff_40, sff_42, sfg1_48, sgf_70, sgf_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_8 * sff_40[k]
                   + f_3 * pc_y[k] * sgf_70[k];

        t_107[k] = f_7 * sff_30[k]
                   + f_3 * pc_z[k] * sgf_70[k];

        t_108[k] = pb_z[k] * sfg0_48[k]
                   - f_6 * pc_z[k] * sfg1_48[k];

        t_109[k] = f_8 * sff_42[k]
                   + f_3 * pc_y[k] * sgf_72[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pb_x, pc_x, sfg0_110, sff_75, sff_76, \
                         sff_77, sff_78, sfg1_110, sgf_76, sgf_77, \
                         sgf_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pb_x[k] * sfg0_110[k]
                   + f_8 * sff_75[k]
                   - f_6 * pc_x[k] * sfg1_110[k];

        t_111[k] = f_7 * sff_76[k]
                   + f_3 * pc_x[k] * sgf_76[k];

        t_112[k] = f_7 * sff_77[k]
                   + f_3 * pc_x[k] * sgf_77[k];

        t_113[k] = f_7 * sff_78[k]
                   + f_3 * pc_x[k] * sgf_78[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pb_x, pc_x, pc_z, sfg0_115, sfg0_117, \
                         sff_36, sff_79, sfg1_115, sfg1_117, sgf_76, \
                         sgf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_7 * sff_79[k]
                   + f_3 * pc_x[k] * sgf_79[k];

        t_115[k] = pb_x[k] * sfg0_115[k]
                   - f_6 * pc_x[k] * sfg1_115[k];

        t_116[k] = f_7 * sff_36[k]
                   + f_3 * pc_z[k] * sgf_76[k];

        t_117[k] = pb_x[k] * sfg0_117[k]
                   - f_6 * pc_x[k] * sfg1_117[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pb_x, pb_y, pc_x, pc_y, sfg0_75, \
                         sfg0_119, sff_49, sff_50, sfg1_75, sfg1_119, sgf_79, \
                         sgf_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_8 * sff_49[k]
                   + f_3 * pc_y[k] * sgf_79[k];

        t_119[k] = pb_x[k] * sfg0_119[k]
                   - f_6 * pc_x[k] * sfg1_119[k];

        t_120[k] = pb_y[k] * sfg0_75[k]
                   - f_6 * pc_y[k] * sfg1_75[k];

        t_121[k] = f_7 * sff_50[k]
                   + f_3 * pc_y[k] * sgf_80[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, pb_x, pc_x, pc_y, pc_z, sfg0_123, sff_40, \
                         sff_52, sff_83, sfg1_123, sgf_80, sgf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_8 * sff_40[k]
                   + f_3 * pc_z[k] * sgf_80[k];

        t_123[k] = pb_x[k] * sfg0_123[k]
                   + f_8 * sff_83[k]
                   - f_6 * pc_x[k] * sfg1_123[k];

        t_124[k] = f_7 * sff_52[k]
                   + f_3 * pc_y[k] * sgf_82[k];
    }
}

static auto
compute_prim_sgg_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sfg0,
                                                          const size_t sff, const size_t sfg1,
                                                          const size_t sgd0, const size_t sgd1,
                                                          const size_t sgf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 1.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfg0_80 = buffer.data(sfg0 + 80);
    const auto *sfg0_90 = buffer.data(sfg0 + 90);
    const auto *sfg0_93 = buffer.data(sfg0 + 93);
    const auto *sfg0_100 = buffer.data(sfg0 + 100);
    const auto *sfg0_102 = buffer.data(sfg0 + 102);
    const auto *sfg0_130 = buffer.data(sfg0 + 130);
    const auto *sfg0_132 = buffer.data(sfg0 + 132);
    const auto *sfg0_134 = buffer.data(sfg0 + 134);
    const auto *sfg0_135 = buffer.data(sfg0 + 135);
    const auto *sfg0_138 = buffer.data(sfg0 + 138);
    const auto *sfg0_140 = buffer.data(sfg0 + 140);
    const auto *sfg0_145 = buffer.data(sfg0 + 145);
    const auto *sfg0_147 = buffer.data(sfg0 + 147);
    const auto *sfg0_149 = buffer.data(sfg0 + 149);

    const auto *sff_46 = buffer.data(sff + 46);
    const auto *sff_50 = buffer.data(sff + 50);
    const auto *sff_56 = buffer.data(sff + 56);
    const auto *sff_59 = buffer.data(sff + 59);
    const auto *sff_60 = buffer.data(sff + 60);
    const auto *sff_62 = buffer.data(sff + 62);
    const auto *sff_66 = buffer.data(sff + 66);
    const auto *sff_67 = buffer.data(sff + 67);
    const auto *sff_68 = buffer.data(sff + 68);
    const auto *sff_69 = buffer.data(sff + 69);
    const auto *sff_70 = buffer.data(sff + 70);
    const auto *sff_72 = buffer.data(sff + 72);
    const auto *sff_76 = buffer.data(sff + 76);
    const auto *sff_79 = buffer.data(sff + 79);
    const auto *sff_80 = buffer.data(sff + 80);
    const auto *sff_82 = buffer.data(sff + 82);
    const auto *sff_86 = buffer.data(sff + 86);
    const auto *sff_87 = buffer.data(sff + 87);
    const auto *sff_88 = buffer.data(sff + 88);
    const auto *sff_89 = buffer.data(sff + 89);
    const auto *sff_90 = buffer.data(sff + 90);
    const auto *sff_92 = buffer.data(sff + 92);
    const auto *sff_93 = buffer.data(sff + 93);
    const auto *sff_95 = buffer.data(sff + 95);
    const auto *sff_96 = buffer.data(sff + 96);
    const auto *sff_97 = buffer.data(sff + 97);
    const auto *sff_98 = buffer.data(sff + 98);
    const auto *sff_99 = buffer.data(sff + 99);

    const auto *sfg1_80 = buffer.data(sfg1 + 80);
    const auto *sfg1_90 = buffer.data(sfg1 + 90);
    const auto *sfg1_93 = buffer.data(sfg1 + 93);
    const auto *sfg1_100 = buffer.data(sfg1 + 100);
    const auto *sfg1_102 = buffer.data(sfg1 + 102);
    const auto *sfg1_130 = buffer.data(sfg1 + 130);
    const auto *sfg1_132 = buffer.data(sfg1 + 132);
    const auto *sfg1_134 = buffer.data(sfg1 + 134);
    const auto *sfg1_135 = buffer.data(sfg1 + 135);
    const auto *sfg1_138 = buffer.data(sfg1 + 138);
    const auto *sfg1_140 = buffer.data(sfg1 + 140);
    const auto *sfg1_145 = buffer.data(sfg1 + 145);
    const auto *sfg1_147 = buffer.data(sfg1 + 147);
    const auto *sfg1_149 = buffer.data(sfg1 + 149);

    const auto *sgd0_60 = buffer.data(sgd0 + 60);
    const auto *sgd0_63 = buffer.data(sgd0 + 63);
    const auto *sgd0_65 = buffer.data(sgd0 + 65);
    const auto *sgd0_71 = buffer.data(sgd0 + 71);
    const auto *sgd0_72 = buffer.data(sgd0 + 72);
    const auto *sgd0_75 = buffer.data(sgd0 + 75);
    const auto *sgd0_77 = buffer.data(sgd0 + 77);
    const auto *sgd0_81 = buffer.data(sgd0 + 81);
    const auto *sgd0_84 = buffer.data(sgd0 + 84);
    const auto *sgd0_87 = buffer.data(sgd0 + 87);
    const auto *sgd0_89 = buffer.data(sgd0 + 89);

    const auto *sgd1_60 = buffer.data(sgd1 + 60);
    const auto *sgd1_63 = buffer.data(sgd1 + 63);
    const auto *sgd1_65 = buffer.data(sgd1 + 65);
    const auto *sgd1_71 = buffer.data(sgd1 + 71);
    const auto *sgd1_72 = buffer.data(sgd1 + 72);
    const auto *sgd1_75 = buffer.data(sgd1 + 75);
    const auto *sgd1_77 = buffer.data(sgd1 + 77);
    const auto *sgd1_81 = buffer.data(sgd1 + 81);
    const auto *sgd1_84 = buffer.data(sgd1 + 84);
    const auto *sgd1_87 = buffer.data(sgd1 + 87);
    const auto *sgd1_89 = buffer.data(sgd1 + 89);

    const auto *sgf_86 = buffer.data(sgf + 86);
    const auto *sgf_87 = buffer.data(sgf + 87);
    const auto *sgf_88 = buffer.data(sgf + 88);
    const auto *sgf_89 = buffer.data(sgf + 89);
    const auto *sgf_90 = buffer.data(sgf + 90);
    const auto *sgf_92 = buffer.data(sgf + 92);
    const auto *sgf_96 = buffer.data(sgf + 96);
    const auto *sgf_97 = buffer.data(sgf + 97);
    const auto *sgf_98 = buffer.data(sgf + 98);
    const auto *sgf_99 = buffer.data(sgf + 99);
    const auto *sgf_100 = buffer.data(sgf + 100);
    const auto *sgf_102 = buffer.data(sgf + 102);
    const auto *sgf_103 = buffer.data(sgf + 103);
    const auto *sgf_105 = buffer.data(sgf + 105);
    const auto *sgf_106 = buffer.data(sgf + 106);
    const auto *sgf_107 = buffer.data(sgf + 107);
    const auto *sgf_108 = buffer.data(sgf + 108);
    const auto *sgf_109 = buffer.data(sgf + 109);
    const auto *sgf_110 = buffer.data(sgf + 110);
    const auto *sgf_112 = buffer.data(sgf + 112);
    const auto *sgf_115 = buffer.data(sgf + 115);
    const auto *sgf_116 = buffer.data(sgf + 116);
    const auto *sgf_117 = buffer.data(sgf + 117);
    const auto *sgf_118 = buffer.data(sgf + 118);
    const auto *sgf_119 = buffer.data(sgf + 119);
    const auto *sgf_120 = buffer.data(sgf + 120);
    const auto *sgf_122 = buffer.data(sgf + 122);
    const auto *sgf_123 = buffer.data(sgf + 123);
    const auto *sgf_125 = buffer.data(sgf + 125);
    const auto *sgf_126 = buffer.data(sgf + 126);
    const auto *sgf_127 = buffer.data(sgf + 127);
    const auto *sgf_128 = buffer.data(sgf + 128);
    const auto *sgf_129 = buffer.data(sgf + 129);
    const auto *sgf_130 = buffer.data(sgf + 130);
    const auto *sgf_132 = buffer.data(sgf + 132);
    const auto *sgf_133 = buffer.data(sgf + 133);
    const auto *sgf_136 = buffer.data(sgf + 136);
    const auto *sgf_137 = buffer.data(sgf + 137);
    const auto *sgf_138 = buffer.data(sgf + 138);
    const auto *sgf_139 = buffer.data(sgf + 139);
    const auto *sgf_140 = buffer.data(sgf + 140);
    const auto *sgf_142 = buffer.data(sgf + 142);
    const auto *sgf_143 = buffer.data(sgf + 143);
    const auto *sgf_145 = buffer.data(sgf + 145);
    const auto *sgf_146 = buffer.data(sgf + 146);
    const auto *sgf_147 = buffer.data(sgf + 147);
    const auto *sgf_148 = buffer.data(sgf + 148);
    const auto *sgf_149 = buffer.data(sgf + 149);

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pb_y, pc_x, pc_y, sfg0_80, sff_86, \
                         sff_87, sff_88, sfg1_80, sgf_86, sgf_87, \
                         sgf_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = pb_y[k] * sfg0_80[k]
                   - f_6 * pc_y[k] * sfg1_80[k];

        t_126[k] = f_7 * sff_86[k]
                   + f_3 * pc_x[k] * sgf_86[k];

        t_127[k] = f_7 * sff_87[k]
                   + f_3 * pc_x[k] * sgf_87[k];

        t_128[k] = f_7 * sff_88[k]
                   + f_3 * pc_x[k] * sgf_88[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, pb_x, pc_x, pc_z, sfg0_130, sfg0_132, \
                         sff_46, sff_89, sfg1_130, sfg1_132, sgf_86, \
                         sgf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_7 * sff_89[k]
                   + f_3 * pc_x[k] * sgf_89[k];

        t_130[k] = pb_x[k] * sfg0_130[k]
                   - f_6 * pc_x[k] * sfg1_130[k];

        t_131[k] = f_8 * sff_46[k]
                   + f_3 * pc_z[k] * sgf_86[k];

        t_132[k] = pb_x[k] * sfg0_132[k]
                   - f_6 * pc_x[k] * sfg1_132[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pb_x, pc_x, pc_y, sfg0_134, sfg0_135, \
                         sff_59, sff_90, sfg1_134, sfg1_135, sgf_89, \
                         sgf_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_7 * sff_59[k]
                   + f_3 * pc_y[k] * sgf_89[k];

        t_134[k] = pb_x[k] * sfg0_134[k]
                   - f_6 * pc_x[k] * sfg1_134[k];

        t_135[k] = pb_x[k] * sfg0_135[k]
                   + f_0 * sff_90[k]
                   - f_6 * pc_x[k] * sfg1_135[k];

        t_136[k] = f_3 * pc_y[k] * sgf_90[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pb_x, pc_x, pc_y, pc_z, sfg0_138, sff_50, \
                         sff_93, sfg1_138, sgf_90, sgf_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_9 * sff_50[k]
                   + f_3 * pc_z[k] * sgf_90[k];

        t_138[k] = pb_x[k] * sfg0_138[k]
                   + f_8 * sff_93[k]
                   - f_6 * pc_x[k] * sfg1_138[k];

        t_139[k] = f_3 * pc_y[k] * sgf_92[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pb_x, pc_x, sfg0_140, sff_95, sff_96, \
                         sff_97, sff_98, sfg1_140, sgf_96, sgf_97, \
                         sgf_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = pb_x[k] * sfg0_140[k]
                   + f_8 * sff_95[k]
                   - f_6 * pc_x[k] * sfg1_140[k];

        t_141[k] = f_7 * sff_96[k]
                   + f_3 * pc_x[k] * sgf_96[k];

        t_142[k] = f_7 * sff_97[k]
                   + f_3 * pc_x[k] * sgf_97[k];

        t_143[k] = f_7 * sff_98[k]
                   + f_3 * pc_x[k] * sgf_98[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pb_x, pc_x, pc_z, sfg0_145, sfg0_147, \
                         sff_56, sff_99, sfg1_145, sfg1_147, sgf_96, \
                         sgf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_7 * sff_99[k]
                   + f_3 * pc_x[k] * sgf_99[k];

        t_145[k] = pb_x[k] * sfg0_145[k]
                   - f_6 * pc_x[k] * sfg1_145[k];

        t_146[k] = f_9 * sff_56[k]
                   + f_3 * pc_z[k] * sgf_96[k];

        t_147[k] = pb_x[k] * sfg0_147[k]
                   - f_6 * pc_x[k] * sfg1_147[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, t_152, pb_x, pc_x, pc_y, pc_z, sfg0_149, \
                         sff_60, sfg1_149, sgd0_60, sgd1_60, sgf_99, \
                         sgf_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_3 * pc_y[k] * sgf_99[k];

        t_149[k] = pb_x[k] * sfg0_149[k]
                   - f_6 * pc_x[k] * sfg1_149[k];

        t_150[k] = f_1 * sgd0_60[k]
                   - f_2 * sgd1_60[k]
                   + f_3 * pc_x[k] * sgf_100[k];

        t_151[k] = f_0 * sff_60[k]
                   + f_3 * pc_y[k] * sgf_100[k];

        t_152[k] = f_3 * pc_z[k] * sgf_100[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, pc_x, pc_y, sff_62, sgd0_63, sgd0_65, \
                         sgd1_63, sgd1_65, sgf_102, sgf_103, sgf_105, \
                         sgf_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_4 * sgd0_63[k]
                   - f_5 * sgd1_63[k]
                   + f_3 * pc_x[k] * sgf_103[k];

        t_154[k] = f_0 * sff_62[k]
                   + f_3 * pc_y[k] * sgf_102[k];

        t_155[k] = f_4 * sgd0_65[k]
                   - f_5 * sgd1_65[k]
                   + f_3 * pc_x[k] * sgf_105[k];

        t_156[k] = f_3 * pc_x[k] * sgf_106[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, pc_x, pc_y, pc_z, sff_66, sgd0_63, \
                         sgd1_63, sgf_106, sgf_107, sgf_108, sgf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_3 * pc_x[k] * sgf_107[k];

        t_158[k] = f_3 * pc_x[k] * sgf_108[k];

        t_159[k] = f_3 * pc_x[k] * sgf_109[k];

        t_160[k] = f_0 * sff_66[k]
                   + f_1 * sgd0_63[k]
                   - f_2 * sgd1_63[k]
                   + f_3 * pc_y[k] * sgf_106[k];

        t_161[k] = f_3 * pc_z[k] * sgf_106[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pb_z, pc_y, pc_z, sfg0_90, sff_68, \
                         sff_69, sfg1_90, sgd0_65, sgd1_65, sgf_108, \
                         sgf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_0 * sff_68[k]
                   + f_4 * sgd0_65[k]
                   - f_5 * sgd1_65[k]
                   + f_3 * pc_y[k] * sgf_108[k];

        t_163[k] = f_0 * sff_69[k]
                   + f_3 * pc_y[k] * sgf_109[k];

        t_164[k] = f_1 * sgd0_65[k]
                   - f_2 * sgd1_65[k]
                   + f_3 * pc_z[k] * sgf_109[k];

        t_165[k] = pb_z[k] * sfg0_90[k]
                   - f_6 * pc_z[k] * sfg1_90[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pb_z, pc_y, pc_z, sfg0_93, sff_60, \
                         sff_70, sff_72, sfg1_93, sgf_110, sgf_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_9 * sff_70[k]
                   + f_3 * pc_y[k] * sgf_110[k];

        t_167[k] = f_7 * sff_60[k]
                   + f_3 * pc_z[k] * sgf_110[k];

        t_168[k] = pb_z[k] * sfg0_93[k]
                   - f_6 * pc_z[k] * sfg1_93[k];

        t_169[k] = f_9 * sff_72[k]
                   + f_3 * pc_y[k] * sgf_112[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, pc_x, sgd0_71, sgd1_71, sgf_115, \
                         sgf_116, sgf_117, sgf_118, sgf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_4 * sgd0_71[k]
                   - f_5 * sgd1_71[k]
                   + f_3 * pc_x[k] * sgf_115[k];

        t_171[k] = f_3 * pc_x[k] * sgf_116[k];

        t_172[k] = f_3 * pc_x[k] * sgf_117[k];

        t_173[k] = f_3 * pc_x[k] * sgf_118[k];

        t_174[k] = f_3 * pc_x[k] * sgf_119[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, pb_z, pc_y, pc_z, sfg0_100, sfg0_102, \
                         sff_66, sff_67, sff_79, sfg1_100, sfg1_102, sgf_116, \
                         sgf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = pb_z[k] * sfg0_100[k]
                   - f_6 * pc_z[k] * sfg1_100[k];

        t_176[k] = f_7 * sff_66[k]
                   + f_3 * pc_z[k] * sgf_116[k];

        t_177[k] = pb_z[k] * sfg0_102[k]
                   + f_8 * sff_67[k]
                   - f_6 * pc_z[k] * sfg1_102[k];

        t_178[k] = f_9 * sff_79[k]
                   + f_3 * pc_y[k] * sgf_119[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, t_182, pc_x, pc_y, pc_z, sff_69, sff_70, sff_80, \
                         sgd0_71, sgd0_72, sgd1_71, sgd1_72, sgf_119, \
                         sgf_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_7 * sff_69[k]
                   + f_1 * sgd0_71[k]
                   - f_2 * sgd1_71[k]
                   + f_3 * pc_z[k] * sgf_119[k];

        t_180[k] = f_1 * sgd0_72[k]
                   - f_2 * sgd1_72[k]
                   + f_3 * pc_x[k] * sgf_120[k];

        t_181[k] = f_8 * sff_80[k]
                   + f_3 * pc_y[k] * sgf_120[k];

        t_182[k] = f_8 * sff_70[k]
                   + f_3 * pc_z[k] * sgf_120[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, pc_x, pc_y, sff_82, sgd0_75, sgd0_77, \
                         sgd1_75, sgd1_77, sgf_122, sgf_123, sgf_125, \
                         sgf_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_4 * sgd0_75[k]
                   - f_5 * sgd1_75[k]
                   + f_3 * pc_x[k] * sgf_123[k];

        t_184[k] = f_8 * sff_82[k]
                   + f_3 * pc_y[k] * sgf_122[k];

        t_185[k] = f_4 * sgd0_77[k]
                   - f_5 * sgd1_77[k]
                   + f_3 * pc_x[k] * sgf_125[k];

        t_186[k] = f_3 * pc_x[k] * sgf_126[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, t_190, t_191, pc_x, pc_y, pc_z, sff_76, sff_86, \
                         sgd0_75, sgd1_75, sgf_126, sgf_127, sgf_128, \
                         sgf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = f_3 * pc_x[k] * sgf_127[k];

        t_188[k] = f_3 * pc_x[k] * sgf_128[k];

        t_189[k] = f_3 * pc_x[k] * sgf_129[k];

        t_190[k] = f_8 * sff_86[k]
                   + f_1 * sgd0_75[k]
                   - f_2 * sgd1_75[k]
                   + f_3 * pc_y[k] * sgf_126[k];

        t_191[k] = f_8 * sff_76[k]
                   + f_3 * pc_z[k] * sgf_126[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pb_y, pc_y, pc_z, sfg0_135, sff_79, \
                         sff_88, sff_89, sfg1_135, sgd0_77, sgd1_77, sgf_128, \
                         sgf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_8 * sff_88[k]
                   + f_4 * sgd0_77[k]
                   - f_5 * sgd1_77[k]
                   + f_3 * pc_y[k] * sgf_128[k];

        t_193[k] = f_8 * sff_89[k]
                   + f_3 * pc_y[k] * sgf_129[k];

        t_194[k] = f_8 * sff_79[k]
                   + f_1 * sgd0_77[k]
                   - f_2 * sgd1_77[k]
                   + f_3 * pc_z[k] * sgf_129[k];

        t_195[k] = pb_y[k] * sfg0_135[k]
                   - f_6 * pc_y[k] * sfg1_135[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pc_x, pc_y, pc_z, sff_80, sff_90, sff_92, \
                         sgd0_81, sgd1_81, sgf_130, sgf_132, sgf_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_7 * sff_90[k]
                   + f_3 * pc_y[k] * sgf_130[k];

        t_197[k] = f_9 * sff_80[k]
                   + f_3 * pc_z[k] * sgf_130[k];

        t_198[k] = f_4 * sgd0_81[k]
                   - f_5 * sgd1_81[k]
                   + f_3 * pc_x[k] * sgf_133[k];

        t_199[k] = f_7 * sff_92[k]
                   + f_3 * pc_y[k] * sgf_132[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, pb_y, pc_x, pc_y, sfg0_140, \
                         sfg1_140, sgf_136, sgf_137, sgf_138, sgf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = pb_y[k] * sfg0_140[k]
                   - f_6 * pc_y[k] * sfg1_140[k];

        t_201[k] = f_3 * pc_x[k] * sgf_136[k];

        t_202[k] = f_3 * pc_x[k] * sgf_137[k];

        t_203[k] = f_3 * pc_x[k] * sgf_138[k];

        t_204[k] = f_3 * pc_x[k] * sgf_139[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, pb_y, pc_y, pc_z, sfg0_145, sfg0_147, sff_86, \
                         sff_96, sff_98, sfg1_145, sfg1_147, sgf_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = pb_y[k] * sfg0_145[k]
                   + f_0 * sff_96[k]
                   - f_6 * pc_y[k] * sfg1_145[k];

        t_206[k] = f_9 * sff_86[k]
                   + f_3 * pc_z[k] * sgf_136[k];

        t_207[k] = pb_y[k] * sfg0_147[k]
                   + f_8 * sff_98[k]
                   - f_6 * pc_y[k] * sfg1_147[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pb_y, pc_x, pc_y, sfg0_149, sff_99, \
                         sfg1_149, sgd0_84, sgd1_84, sgf_139, sgf_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_7 * sff_99[k]
                   + f_3 * pc_y[k] * sgf_139[k];

        t_209[k] = pb_y[k] * sfg0_149[k]
                   - f_6 * pc_y[k] * sfg1_149[k];

        t_210[k] = f_1 * sgd0_84[k]
                   - f_2 * sgd1_84[k]
                   + f_3 * pc_x[k] * sgf_140[k];

        t_211[k] = f_3 * pc_y[k] * sgf_140[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pc_x, pc_y, pc_z, sff_90, sgd0_87, \
                         sgd0_89, sgd1_87, sgd1_89, sgf_140, sgf_142, sgf_143, \
                         sgf_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_0 * sff_90[k]
                   + f_3 * pc_z[k] * sgf_140[k];

        t_213[k] = f_4 * sgd0_87[k]
                   - f_5 * sgd1_87[k]
                   + f_3 * pc_x[k] * sgf_143[k];

        t_214[k] = f_3 * pc_y[k] * sgf_142[k];

        t_215[k] = f_4 * sgd0_89[k]
                   - f_5 * sgd1_89[k]
                   + f_3 * pc_x[k] * sgf_145[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, t_219, t_220, t_221, pc_x, pc_y, pc_z, sff_96, \
                         sgd0_87, sgd1_87, sgf_146, sgf_147, sgf_148, \
                         sgf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_3 * pc_x[k] * sgf_146[k];

        t_217[k] = f_3 * pc_x[k] * sgf_147[k];

        t_218[k] = f_3 * pc_x[k] * sgf_148[k];

        t_219[k] = f_3 * pc_x[k] * sgf_149[k];

        t_220[k] = f_1 * sgd0_87[k]
                   - f_2 * sgd1_87[k]
                   + f_3 * pc_y[k] * sgf_146[k];

        t_221[k] = f_0 * sff_96[k]
                   + f_3 * pc_z[k] * sgf_146[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, pc_y, pc_z, sff_99, sgd0_89, sgd1_89, sgf_148, \
                         sgf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_4 * sgd0_89[k]
                   - f_5 * sgd1_89[k]
                   + f_3 * pc_y[k] * sgf_148[k];

        t_223[k] = f_3 * pc_y[k] * sgf_149[k];

        t_224[k] = f_0 * sff_99[k]
                   + f_1 * sgd0_89[k]
                   - f_2 * sgd1_89[k]
                   + f_3 * pc_z[k] * sgf_149[k];
    }
}

auto
compute_prim_sgg_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sfg0, const size_t sff,
                                                   const size_t sfg1, const size_t sgd0,
                                                   const size_t sgd1, const size_t sgf,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_sgg_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, sfg0, sff,
                                                              sfg1, sgd0, sgd1, sgf, ncols,
                                                              gamma, p, q);

    compute_prim_sgg_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, sfg0, sff,
                                                              sfg1, sgd0, sgd1, sgf, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
