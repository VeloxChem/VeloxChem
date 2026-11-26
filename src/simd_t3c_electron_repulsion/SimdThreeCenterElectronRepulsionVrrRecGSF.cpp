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


#include "SimdThreeCenterElectronRepulsionVrrRecGSF.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_gsf_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t fsf0,
                                                          const size_t fsd, const size_t fsf1,
                                                          const size_t gsp0, const size_t gsp1,
                                                          const size_t gsd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 1.5 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_9 = 1.0 / q;

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
    auto *t_129 = buffer.data(target + 129);
    auto *t_130 = buffer.data(target + 130);
    auto *t_131 = buffer.data(target + 131);
    auto *t_132 = buffer.data(target + 132);
    auto *t_133 = buffer.data(target + 133);
    auto *t_134 = buffer.data(target + 134);
    auto *t_135 = buffer.data(target + 135);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *fsf0_0 = buffer.data(fsf0 + 0);
    const auto *fsf0_6 = buffer.data(fsf0 + 6);
    const auto *fsf0_9 = buffer.data(fsf0 + 9);
    const auto *fsf0_16 = buffer.data(fsf0 + 16);
    const auto *fsf0_20 = buffer.data(fsf0 + 20);
    const auto *fsf0_29 = buffer.data(fsf0 + 29);
    const auto *fsf0_30 = buffer.data(fsf0 + 30);
    const auto *fsf0_50 = buffer.data(fsf0 + 50);
    const auto *fsf0_60 = buffer.data(fsf0 + 60);
    const auto *fsf0_61 = buffer.data(fsf0 + 61);
    const auto *fsf0_66 = buffer.data(fsf0 + 66);
    const auto *fsf0_69 = buffer.data(fsf0 + 69);
    const auto *fsf0_76 = buffer.data(fsf0 + 76);
    const auto *fsf0_79 = buffer.data(fsf0 + 79);
    const auto *fsf0_86 = buffer.data(fsf0 + 86);
    const auto *fsf0_89 = buffer.data(fsf0 + 89);
    const auto *fsf0_90 = buffer.data(fsf0 + 90);
    const auto *fsf0_92 = buffer.data(fsf0 + 92);
    const auto *fsf0_96 = buffer.data(fsf0 + 96);
    const auto *fsf0_97 = buffer.data(fsf0 + 97);
    const auto *fsf0_99 = buffer.data(fsf0 + 99);

    const auto *fsd_0 = buffer.data(fsd + 0);
    const auto *fsd_3 = buffer.data(fsd + 3);
    const auto *fsd_5 = buffer.data(fsd + 5);
    const auto *fsd_6 = buffer.data(fsd + 6);
    const auto *fsd_9 = buffer.data(fsd + 9);
    const auto *fsd_11 = buffer.data(fsd + 11);
    const auto *fsd_12 = buffer.data(fsd + 12);
    const auto *fsd_15 = buffer.data(fsd + 15);
    const auto *fsd_17 = buffer.data(fsd + 17);
    const auto *fsd_18 = buffer.data(fsd + 18);
    const auto *fsd_21 = buffer.data(fsd + 21);
    const auto *fsd_23 = buffer.data(fsd + 23);
    const auto *fsd_24 = buffer.data(fsd + 24);
    const auto *fsd_27 = buffer.data(fsd + 27);
    const auto *fsd_28 = buffer.data(fsd + 28);
    const auto *fsd_29 = buffer.data(fsd + 29);
    const auto *fsd_30 = buffer.data(fsd + 30);
    const auto *fsd_33 = buffer.data(fsd + 33);
    const auto *fsd_35 = buffer.data(fsd + 35);
    const auto *fsd_36 = buffer.data(fsd + 36);
    const auto *fsd_39 = buffer.data(fsd + 39);
    const auto *fsd_41 = buffer.data(fsd + 41);
    const auto *fsd_45 = buffer.data(fsd + 45);
    const auto *fsd_46 = buffer.data(fsd + 46);
    const auto *fsd_47 = buffer.data(fsd + 47);
    const auto *fsd_51 = buffer.data(fsd + 51);
    const auto *fsd_52 = buffer.data(fsd + 52);
    const auto *fsd_53 = buffer.data(fsd + 53);
    const auto *fsd_54 = buffer.data(fsd + 54);
    const auto *fsd_57 = buffer.data(fsd + 57);
    const auto *fsd_59 = buffer.data(fsd + 59);

    const auto *fsf1_0 = buffer.data(fsf1 + 0);
    const auto *fsf1_6 = buffer.data(fsf1 + 6);
    const auto *fsf1_9 = buffer.data(fsf1 + 9);
    const auto *fsf1_16 = buffer.data(fsf1 + 16);
    const auto *fsf1_20 = buffer.data(fsf1 + 20);
    const auto *fsf1_29 = buffer.data(fsf1 + 29);
    const auto *fsf1_30 = buffer.data(fsf1 + 30);
    const auto *fsf1_50 = buffer.data(fsf1 + 50);
    const auto *fsf1_60 = buffer.data(fsf1 + 60);
    const auto *fsf1_61 = buffer.data(fsf1 + 61);
    const auto *fsf1_66 = buffer.data(fsf1 + 66);
    const auto *fsf1_69 = buffer.data(fsf1 + 69);
    const auto *fsf1_76 = buffer.data(fsf1 + 76);
    const auto *fsf1_79 = buffer.data(fsf1 + 79);
    const auto *fsf1_86 = buffer.data(fsf1 + 86);
    const auto *fsf1_89 = buffer.data(fsf1 + 89);
    const auto *fsf1_90 = buffer.data(fsf1 + 90);
    const auto *fsf1_92 = buffer.data(fsf1 + 92);
    const auto *fsf1_96 = buffer.data(fsf1 + 96);
    const auto *fsf1_97 = buffer.data(fsf1 + 97);
    const auto *fsf1_99 = buffer.data(fsf1 + 99);

    const auto *gsp0_0 = buffer.data(gsp0 + 0);
    const auto *gsp0_1 = buffer.data(gsp0 + 1);
    const auto *gsp0_2 = buffer.data(gsp0 + 2);
    const auto *gsp0_4 = buffer.data(gsp0 + 4);
    const auto *gsp0_8 = buffer.data(gsp0 + 8);
    const auto *gsp0_9 = buffer.data(gsp0 + 9);
    const auto *gsp0_10 = buffer.data(gsp0 + 10);
    const auto *gsp0_11 = buffer.data(gsp0 + 11);
    const auto *gsp0_15 = buffer.data(gsp0 + 15);
    const auto *gsp0_16 = buffer.data(gsp0 + 16);
    const auto *gsp0_17 = buffer.data(gsp0 + 17);
    const auto *gsp0_30 = buffer.data(gsp0 + 30);
    const auto *gsp0_31 = buffer.data(gsp0 + 31);
    const auto *gsp0_32 = buffer.data(gsp0 + 32);
    const auto *gsp0_35 = buffer.data(gsp0 + 35);
    const auto *gsp0_36 = buffer.data(gsp0 + 36);
    const auto *gsp0_37 = buffer.data(gsp0 + 37);
    const auto *gsp0_38 = buffer.data(gsp0 + 38);
    const auto *gsp0_40 = buffer.data(gsp0 + 40);

    const auto *gsp1_0 = buffer.data(gsp1 + 0);
    const auto *gsp1_1 = buffer.data(gsp1 + 1);
    const auto *gsp1_2 = buffer.data(gsp1 + 2);
    const auto *gsp1_4 = buffer.data(gsp1 + 4);
    const auto *gsp1_8 = buffer.data(gsp1 + 8);
    const auto *gsp1_9 = buffer.data(gsp1 + 9);
    const auto *gsp1_10 = buffer.data(gsp1 + 10);
    const auto *gsp1_11 = buffer.data(gsp1 + 11);
    const auto *gsp1_15 = buffer.data(gsp1 + 15);
    const auto *gsp1_16 = buffer.data(gsp1 + 16);
    const auto *gsp1_17 = buffer.data(gsp1 + 17);
    const auto *gsp1_30 = buffer.data(gsp1 + 30);
    const auto *gsp1_31 = buffer.data(gsp1 + 31);
    const auto *gsp1_32 = buffer.data(gsp1 + 32);
    const auto *gsp1_35 = buffer.data(gsp1 + 35);
    const auto *gsp1_36 = buffer.data(gsp1 + 36);
    const auto *gsp1_37 = buffer.data(gsp1 + 37);
    const auto *gsp1_38 = buffer.data(gsp1 + 38);
    const auto *gsp1_40 = buffer.data(gsp1 + 40);

    const auto *gsd_0 = buffer.data(gsd + 0);
    const auto *gsd_2 = buffer.data(gsd + 2);
    const auto *gsd_3 = buffer.data(gsd + 3);
    const auto *gsd_5 = buffer.data(gsd + 5);
    const auto *gsd_6 = buffer.data(gsd + 6);
    const auto *gsd_7 = buffer.data(gsd + 7);
    const auto *gsd_9 = buffer.data(gsd + 9);
    const auto *gsd_11 = buffer.data(gsd + 11);
    const auto *gsd_12 = buffer.data(gsd + 12);
    const auto *gsd_14 = buffer.data(gsd + 14);
    const auto *gsd_15 = buffer.data(gsd + 15);
    const auto *gsd_16 = buffer.data(gsd + 16);
    const auto *gsd_17 = buffer.data(gsd + 17);
    const auto *gsd_18 = buffer.data(gsd + 18);
    const auto *gsd_19 = buffer.data(gsd + 19);
    const auto *gsd_21 = buffer.data(gsd + 21);
    const auto *gsd_23 = buffer.data(gsd + 23);
    const auto *gsd_24 = buffer.data(gsd + 24);
    const auto *gsd_27 = buffer.data(gsd + 27);
    const auto *gsd_28 = buffer.data(gsd + 28);
    const auto *gsd_29 = buffer.data(gsd + 29);
    const auto *gsd_30 = buffer.data(gsd + 30);
    const auto *gsd_32 = buffer.data(gsd + 32);
    const auto *gsd_33 = buffer.data(gsd + 33);
    const auto *gsd_34 = buffer.data(gsd + 34);
    const auto *gsd_35 = buffer.data(gsd + 35);
    const auto *gsd_36 = buffer.data(gsd + 36);
    const auto *gsd_37 = buffer.data(gsd + 37);
    const auto *gsd_39 = buffer.data(gsd + 39);
    const auto *gsd_41 = buffer.data(gsd + 41);
    const auto *gsd_42 = buffer.data(gsd + 42);
    const auto *gsd_45 = buffer.data(gsd + 45);
    const auto *gsd_46 = buffer.data(gsd + 46);
    const auto *gsd_47 = buffer.data(gsd + 47);
    const auto *gsd_48 = buffer.data(gsd + 48);
    const auto *gsd_51 = buffer.data(gsd + 51);
    const auto *gsd_52 = buffer.data(gsd + 52);
    const auto *gsd_53 = buffer.data(gsd + 53);
    const auto *gsd_54 = buffer.data(gsd + 54);
    const auto *gsd_56 = buffer.data(gsd + 56);
    const auto *gsd_57 = buffer.data(gsd + 57);
    const auto *gsd_59 = buffer.data(gsd + 59);
    const auto *gsd_60 = buffer.data(gsd + 60);
    const auto *gsd_61 = buffer.data(gsd + 61);
    const auto *gsd_63 = buffer.data(gsd + 63);
    const auto *gsd_64 = buffer.data(gsd + 64);
    const auto *gsd_65 = buffer.data(gsd + 65);
    const auto *gsd_68 = buffer.data(gsd + 68);
    const auto *gsd_69 = buffer.data(gsd + 69);
    const auto *gsd_70 = buffer.data(gsd + 70);
    const auto *gsd_71 = buffer.data(gsd + 71);
    const auto *gsd_72 = buffer.data(gsd + 72);
    const auto *gsd_73 = buffer.data(gsd + 73);
    const auto *gsd_74 = buffer.data(gsd + 74);
    const auto *gsd_75 = buffer.data(gsd + 75);
    const auto *gsd_76 = buffer.data(gsd + 76);
    const auto *gsd_77 = buffer.data(gsd + 77);
    const auto *gsd_79 = buffer.data(gsd + 79);
    const auto *gsd_81 = buffer.data(gsd + 81);
    const auto *gsd_82 = buffer.data(gsd + 82);
    const auto *gsd_83 = buffer.data(gsd + 83);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, pc_z, fsd_0, fsd_3, gsp0_0, \
                         gsp1_0, gsd_0, gsd_2, gsd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fsd_0[k]
                 + f_1 * gsp0_0[k]
                 - f_2 * gsp1_0[k]
                 + f_3 * pc_x[k] * gsd_0[k];

        t_1[k] = f_3 * pc_y[k] * gsd_0[k];

        t_2[k] = f_3 * pc_z[k] * gsd_0[k];

        t_3[k] = f_0 * fsd_3[k]
                 + f_3 * pc_x[k] * gsd_3[k];

        t_4[k] = f_3 * pc_y[k] * gsd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pc_x, pc_y, pc_z, fsd_5, gsp0_1, gsp0_2, \
                         gsp1_1, gsp1_2, gsd_3, gsd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * fsd_5[k]
                 + f_3 * pc_x[k] * gsd_5[k];

        t_6[k] = f_1 * gsp0_1[k]
                 - f_2 * gsp1_1[k]
                 + f_3 * pc_y[k] * gsd_3[k];

        t_7[k] = f_3 * pc_z[k] * gsd_3[k];

        t_8[k] = f_3 * pc_y[k] * gsd_5[k];

        t_9[k] = f_1 * gsp0_2[k]
                 - f_2 * gsp1_2[k]
                 + f_3 * pc_z[k] * gsd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pc_x, pc_y, pc_z, fsf0_0, fsd_0, \
                         fsd_9, fsf1_0, gsd_6, gsd_7, gsd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * fsf0_0[k]
                  - f_4 * pc_y[k] * fsf1_0[k];

        t_11[k] = f_5 * fsd_0[k]
                  + f_3 * pc_y[k] * gsd_6[k];

        t_12[k] = f_3 * pc_z[k] * gsd_6[k];

        t_13[k] = f_6 * fsd_9[k]
                  + f_3 * pc_x[k] * gsd_9[k];

        t_14[k] = f_3 * pc_z[k] * gsd_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pc_x, pc_y, pc_z, fsd_3, fsd_5, fsd_11, \
                         gsp0_4, gsp1_4, gsd_9, gsd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_6 * fsd_11[k]
                  + f_3 * pc_x[k] * gsd_11[k];

        t_16[k] = f_5 * fsd_3[k]
                  + f_1 * gsp0_4[k]
                  - f_2 * gsp1_4[k]
                  + f_3 * pc_y[k] * gsd_9[k];

        t_17[k] = f_3 * pc_z[k] * gsd_9[k];

        t_18[k] = f_5 * fsd_5[k]
                  + f_3 * pc_y[k] * gsd_11[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_y, pa_z, pc_y, pc_z, fsf0_0, fsf0_9, \
                         fsd_0, fsf1_0, fsf1_9, gsd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_y[k] * fsf0_9[k]
                  - f_4 * pc_y[k] * fsf1_9[k];

        t_20[k] = pa_z[k] * fsf0_0[k]
                  - f_4 * pc_z[k] * fsf1_0[k];

        t_21[k] = f_3 * pc_y[k] * gsd_12[k];

        t_22[k] = f_5 * fsd_0[k]
                  + f_3 * pc_z[k] * gsd_12[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_z, pc_x, pc_y, pc_z, fsf0_6, fsd_15, \
                         fsd_17, fsf1_6, gsd_14, gsd_15, gsd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_6 * fsd_15[k]
                  + f_3 * pc_x[k] * gsd_15[k];

        t_24[k] = f_3 * pc_y[k] * gsd_14[k];

        t_25[k] = f_6 * fsd_17[k]
                  + f_3 * pc_x[k] * gsd_17[k];

        t_26[k] = pa_z[k] * fsf0_6[k]
                  - f_4 * pc_z[k] * fsf1_6[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pc_x, pc_y, pc_z, fsd_5, fsd_18, gsp0_8, \
                         gsp0_9, gsp1_8, gsp1_9, gsd_16, gsd_17, \
                         gsd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_7 * gsp0_8[k]
                  - f_8 * gsp1_8[k]
                  + f_3 * pc_y[k] * gsd_16[k];

        t_28[k] = f_3 * pc_y[k] * gsd_17[k];

        t_29[k] = f_5 * fsd_5[k]
                  + f_1 * gsp0_8[k]
                  - f_2 * gsp1_8[k]
                  + f_3 * pc_z[k] * gsd_17[k];

        t_30[k] = f_9 * fsd_18[k]
                  + f_1 * gsp0_9[k]
                  - f_2 * gsp1_9[k]
                  + f_3 * pc_x[k] * gsd_18[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pc_x, pc_y, pc_z, fsd_6, fsd_21, \
                         fsd_23, gsd_18, gsd_19, gsd_21, gsd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_9 * fsd_6[k]
                  + f_3 * pc_y[k] * gsd_18[k];

        t_32[k] = f_3 * pc_z[k] * gsd_18[k];

        t_33[k] = f_9 * fsd_21[k]
                  + f_3 * pc_x[k] * gsd_21[k];

        t_34[k] = f_3 * pc_z[k] * gsd_19[k];

        t_35[k] = f_9 * fsd_23[k]
                  + f_3 * pc_x[k] * gsd_23[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pc_y, pc_z, fsd_9, fsd_11, gsp0_10, gsp0_11, \
                         gsp1_10, gsp1_11, gsd_21, gsd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_9 * fsd_9[k]
                  + f_1 * gsp0_10[k]
                  - f_2 * gsp1_10[k]
                  + f_3 * pc_y[k] * gsd_21[k];

        t_37[k] = f_3 * pc_z[k] * gsd_21[k];

        t_38[k] = f_9 * fsd_11[k]
                  + f_3 * pc_y[k] * gsd_23[k];

        t_39[k] = f_1 * gsp0_11[k]
                  - f_2 * gsp1_11[k]
                  + f_3 * pc_z[k] * gsd_23[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pc_x, pc_y, pc_z, fsf0_20, fsd_6, \
                         fsd_12, fsd_27, fsf1_20, gsd_24, gsd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_y[k] * fsf0_20[k]
                  - f_4 * pc_y[k] * fsf1_20[k];

        t_41[k] = f_5 * fsd_12[k]
                  + f_3 * pc_y[k] * gsd_24[k];

        t_42[k] = f_5 * fsd_6[k]
                  + f_3 * pc_z[k] * gsd_24[k];

        t_43[k] = f_9 * fsd_27[k]
                  + f_3 * pc_x[k] * gsd_27[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_z, pc_x, pc_z, fsf0_16, fsd_9, fsd_28, \
                         fsd_29, fsf1_16, gsd_27, gsd_28, gsd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_9 * fsd_28[k]
                  + f_3 * pc_x[k] * gsd_28[k];

        t_45[k] = f_9 * fsd_29[k]
                  + f_3 * pc_x[k] * gsd_29[k];

        t_46[k] = pa_z[k] * fsf0_16[k]
                  - f_4 * pc_z[k] * fsf1_16[k];

        t_47[k] = f_5 * fsd_9[k]
                  + f_3 * pc_z[k] * gsd_27[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_y, pc_x, pc_y, fsf0_29, fsd_17, fsd_30, \
                         fsf1_29, gsp0_15, gsp1_15, gsd_29, gsd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_5 * fsd_17[k]
                  + f_3 * pc_y[k] * gsd_29[k];

        t_49[k] = pa_y[k] * fsf0_29[k]
                  - f_4 * pc_y[k] * fsf1_29[k];

        t_50[k] = f_9 * fsd_30[k]
                  + f_1 * gsp0_15[k]
                  - f_2 * gsp1_15[k]
                  + f_3 * pc_x[k] * gsd_30[k];

        t_51[k] = f_3 * pc_y[k] * gsd_30[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pc_x, pc_y, pc_z, fsd_12, fsd_33, fsd_35, \
                         gsd_30, gsd_32, gsd_33, gsd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_9 * fsd_12[k]
                  + f_3 * pc_z[k] * gsd_30[k];

        t_53[k] = f_9 * fsd_33[k]
                  + f_3 * pc_x[k] * gsd_33[k];

        t_54[k] = f_3 * pc_y[k] * gsd_32[k];

        t_55[k] = f_9 * fsd_35[k]
                  + f_3 * pc_x[k] * gsd_35[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pc_y, pc_z, fsd_17, gsp0_16, gsp0_17, \
                         gsp1_16, gsp1_17, gsd_33, gsd_34, gsd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_1 * gsp0_16[k]
                  - f_2 * gsp1_16[k]
                  + f_3 * pc_y[k] * gsd_33[k];

        t_57[k] = f_7 * gsp0_17[k]
                  - f_8 * gsp1_17[k]
                  + f_3 * pc_y[k] * gsd_34[k];

        t_58[k] = f_3 * pc_y[k] * gsd_35[k];

        t_59[k] = f_9 * fsd_17[k]
                  + f_1 * gsp0_17[k]
                  - f_2 * gsp1_17[k]
                  + f_3 * pc_z[k] * gsd_35[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_x, pc_x, pc_y, pc_z, fsf0_60, fsd_18, \
                         fsd_36, fsd_39, fsf1_60, gsd_36, gsd_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pa_x[k] * fsf0_60[k]
                  + f_6 * fsd_36[k]
                  - f_4 * pc_x[k] * fsf1_60[k];

        t_61[k] = f_6 * fsd_18[k]
                  + f_3 * pc_y[k] * gsd_36[k];

        t_62[k] = f_3 * pc_z[k] * gsd_36[k];

        t_63[k] = f_5 * fsd_39[k]
                  + f_3 * pc_x[k] * gsd_39[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, t_68, pa_x, pc_x, pc_y, pc_z, fsf0_66, \
                         fsd_23, fsd_41, fsf1_66, gsd_37, gsd_39, \
                         gsd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_3 * pc_z[k] * gsd_37[k];

        t_65[k] = f_5 * fsd_41[k]
                  + f_3 * pc_x[k] * gsd_41[k];

        t_66[k] = pa_x[k] * fsf0_66[k]
                  - f_4 * pc_x[k] * fsf1_66[k];

        t_67[k] = f_3 * pc_z[k] * gsd_39[k];

        t_68[k] = f_6 * fsd_23[k]
                  + f_3 * pc_y[k] * gsd_41[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_x, pa_z, pc_x, pc_y, pc_z, fsf0_30, \
                         fsf0_69, fsd_18, fsd_24, fsf1_30, fsf1_69, \
                         gsd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = pa_x[k] * fsf0_69[k]
                  - f_4 * pc_x[k] * fsf1_69[k];

        t_70[k] = pa_z[k] * fsf0_30[k]
                  - f_4 * pc_z[k] * fsf1_30[k];

        t_71[k] = f_9 * fsd_24[k]
                  + f_3 * pc_y[k] * gsd_42[k];

        t_72[k] = f_5 * fsd_18[k]
                  + f_3 * pc_z[k] * gsd_42[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_x, pc_x, fsf0_76, fsd_45, fsd_46, fsd_47, \
                         fsf1_76, gsd_45, gsd_46, gsd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_5 * fsd_45[k]
                  + f_3 * pc_x[k] * gsd_45[k];

        t_74[k] = f_5 * fsd_46[k]
                  + f_3 * pc_x[k] * gsd_46[k];

        t_75[k] = f_5 * fsd_47[k]
                  + f_3 * pc_x[k] * gsd_47[k];

        t_76[k] = pa_x[k] * fsf0_76[k]
                  - f_4 * pc_x[k] * fsf1_76[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pa_x, pc_x, pc_y, pc_z, fsf0_79, fsd_21, fsd_29, \
                         fsf1_79, gsd_45, gsd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_5 * fsd_21[k]
                  + f_3 * pc_z[k] * gsd_45[k];

        t_78[k] = f_9 * fsd_29[k]
                  + f_3 * pc_y[k] * gsd_47[k];

        t_79[k] = pa_x[k] * fsf0_79[k]
                  - f_4 * pc_x[k] * fsf1_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pa_y, pc_x, pc_y, pc_z, fsf0_50, fsd_24, \
                         fsd_30, fsd_51, fsf1_50, gsd_48, gsd_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = pa_y[k] * fsf0_50[k]
                  - f_4 * pc_y[k] * fsf1_50[k];

        t_81[k] = f_5 * fsd_30[k]
                  + f_3 * pc_y[k] * gsd_48[k];

        t_82[k] = f_9 * fsd_24[k]
                  + f_3 * pc_z[k] * gsd_48[k];

        t_83[k] = f_5 * fsd_51[k]
                  + f_3 * pc_x[k] * gsd_51[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pa_x, pc_x, pc_z, fsf0_86, fsd_27, fsd_52, \
                         fsd_53, fsf1_86, gsd_51, gsd_52, gsd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_5 * fsd_52[k]
                  + f_3 * pc_x[k] * gsd_52[k];

        t_85[k] = f_5 * fsd_53[k]
                  + f_3 * pc_x[k] * gsd_53[k];

        t_86[k] = pa_x[k] * fsf0_86[k]
                  - f_4 * pc_x[k] * fsf1_86[k];

        t_87[k] = f_9 * fsd_27[k]
                  + f_3 * pc_z[k] * gsd_51[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_x, pc_x, pc_y, fsf0_89, fsf0_90, fsd_35, \
                         fsd_54, fsf1_89, fsf1_90, gsd_53, gsd_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_5 * fsd_35[k]
                  + f_3 * pc_y[k] * gsd_53[k];

        t_89[k] = pa_x[k] * fsf0_89[k]
                  - f_4 * pc_x[k] * fsf1_89[k];

        t_90[k] = pa_x[k] * fsf0_90[k]
                  + f_6 * fsd_54[k]
                  - f_4 * pc_x[k] * fsf1_90[k];

        t_91[k] = f_3 * pc_y[k] * gsd_54[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pc_x, pc_y, pc_z, fsd_30, fsd_57, fsd_59, \
                         gsd_54, gsd_56, gsd_57, gsd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_6 * fsd_30[k]
                  + f_3 * pc_z[k] * gsd_54[k];

        t_93[k] = f_5 * fsd_57[k]
                  + f_3 * pc_x[k] * gsd_57[k];

        t_94[k] = f_3 * pc_y[k] * gsd_56[k];

        t_95[k] = f_5 * fsd_59[k]
                  + f_3 * pc_x[k] * gsd_59[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pa_x, pc_x, pc_y, fsf0_96, fsf0_97, fsf0_99, \
                         fsf1_96, fsf1_97, fsf1_99, gsd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pa_x[k] * fsf0_96[k]
                  - f_4 * pc_x[k] * fsf1_96[k];

        t_97[k] = pa_x[k] * fsf0_97[k]
                  - f_4 * pc_x[k] * fsf1_97[k];

        t_98[k] = f_3 * pc_y[k] * gsd_59[k];

        t_99[k] = pa_x[k] * fsf0_99[k]
                  - f_4 * pc_x[k] * fsf1_99[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, pc_x, pc_z, gsp0_30, gsp0_31, \
                         gsp1_30, gsp1_31, gsd_60, gsd_61, gsd_63, \
                         gsd_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_1 * gsp0_30[k]
                   - f_2 * gsp1_30[k]
                   + f_3 * pc_x[k] * gsd_60[k];

        t_101[k] = f_7 * gsp0_31[k]
                   - f_8 * gsp1_31[k]
                   + f_3 * pc_x[k] * gsd_61[k];

        t_102[k] = f_3 * pc_z[k] * gsd_60[k];

        t_103[k] = f_3 * pc_x[k] * gsd_63[k];

        t_104[k] = f_3 * pc_x[k] * gsd_64[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, t_109, pc_x, pc_y, pc_z, fsd_39, fsd_41, \
                         gsp0_31, gsp0_32, gsp1_31, gsp1_32, gsd_63, \
                         gsd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_3 * pc_x[k] * gsd_65[k];

        t_106[k] = f_0 * fsd_39[k]
                   + f_1 * gsp0_31[k]
                   - f_2 * gsp1_31[k]
                   + f_3 * pc_y[k] * gsd_63[k];

        t_107[k] = f_3 * pc_z[k] * gsd_63[k];

        t_108[k] = f_0 * fsd_41[k]
                   + f_3 * pc_y[k] * gsd_65[k];

        t_109[k] = f_1 * gsp0_32[k]
                   - f_2 * gsp1_32[k]
                   + f_3 * pc_z[k] * gsd_65[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pa_z, pc_x, pc_z, fsf0_60, fsf0_61, \
                         fsf1_60, fsf1_61, gsp0_35, gsp1_35, gsd_68, \
                         gsd_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pa_z[k] * fsf0_60[k]
                   - f_4 * pc_z[k] * fsf1_60[k];

        t_111[k] = pa_z[k] * fsf0_61[k]
                   - f_4 * pc_z[k] * fsf1_61[k];

        t_112[k] = f_7 * gsp0_35[k]
                   - f_8 * gsp1_35[k]
                   + f_3 * pc_x[k] * gsd_68[k];

        t_113[k] = f_3 * pc_x[k] * gsd_69[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, t_118, pa_z, pc_x, pc_y, pc_z, fsf0_66, \
                         fsd_39, fsd_47, fsf1_66, gsd_69, gsd_70, \
                         gsd_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_3 * pc_x[k] * gsd_70[k];

        t_115[k] = f_3 * pc_x[k] * gsd_71[k];

        t_116[k] = pa_z[k] * fsf0_66[k]
                   - f_4 * pc_z[k] * fsf1_66[k];

        t_117[k] = f_5 * fsd_39[k]
                   + f_3 * pc_z[k] * gsd_69[k];

        t_118[k] = f_6 * fsd_47[k]
                   + f_3 * pc_y[k] * gsd_71[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, pc_x, pc_z, fsd_41, gsp0_35, gsp0_36, gsp0_37, \
                         gsp1_35, gsp1_36, gsp1_37, gsd_71, gsd_72, \
                         gsd_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = f_5 * fsd_41[k]
                   + f_1 * gsp0_35[k]
                   - f_2 * gsp1_35[k]
                   + f_3 * pc_z[k] * gsd_71[k];

        t_120[k] = f_1 * gsp0_36[k]
                   - f_2 * gsp1_36[k]
                   + f_3 * pc_x[k] * gsd_72[k];

        t_121[k] = f_7 * gsp0_37[k]
                   - f_8 * gsp1_37[k]
                   + f_3 * pc_x[k] * gsd_73[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, pc_x, pc_y, fsd_51, gsp0_37, \
                         gsp0_38, gsp1_37, gsp1_38, gsd_74, gsd_75, gsd_76, \
                         gsd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_7 * gsp0_38[k]
                   - f_8 * gsp1_38[k]
                   + f_3 * pc_x[k] * gsd_74[k];

        t_123[k] = f_3 * pc_x[k] * gsd_75[k];

        t_124[k] = f_3 * pc_x[k] * gsd_76[k];

        t_125[k] = f_3 * pc_x[k] * gsd_77[k];

        t_126[k] = f_9 * fsd_51[k]
                   + f_1 * gsp0_37[k]
                   - f_2 * gsp1_37[k]
                   + f_3 * pc_y[k] * gsd_75[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, pa_y, pc_y, pc_z, fsf0_90, fsd_45, \
                         fsd_47, fsd_53, fsf1_90, gsp0_38, gsp1_38, gsd_75, \
                         gsd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_9 * fsd_45[k]
                   + f_3 * pc_z[k] * gsd_75[k];

        t_128[k] = f_9 * fsd_53[k]
                   + f_3 * pc_y[k] * gsd_77[k];

        t_129[k] = f_9 * fsd_47[k]
                   + f_1 * gsp0_38[k]
                   - f_2 * gsp1_38[k]
                   + f_3 * pc_z[k] * gsd_77[k];

        t_130[k] = pa_y[k] * fsf0_90[k]
                   - f_4 * pc_y[k] * fsf1_90[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, t_135, pa_y, pc_x, pc_y, fsf0_92, \
                         fsf1_92, gsp0_40, gsp1_40, gsd_79, gsd_81, gsd_82, \
                         gsd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_7 * gsp0_40[k]
                   - f_8 * gsp1_40[k]
                   + f_3 * pc_x[k] * gsd_79[k];

        t_132[k] = pa_y[k] * fsf0_92[k]
                   - f_4 * pc_y[k] * fsf1_92[k];

        t_133[k] = f_3 * pc_x[k] * gsd_81[k];

        t_134[k] = f_3 * pc_x[k] * gsd_82[k];

        t_135[k] = f_3 * pc_x[k] * gsd_83[k];
    }
}

static auto
compute_prim_gsf_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t fsf0,
                                                          const size_t fsd, const size_t fsf1,
                                                          const size_t gsp0, const size_t gsp1,
                                                          const size_t gsd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 1.5 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *fsf0_96 = buffer.data(fsf0 + 96);
    const auto *fsf0_99 = buffer.data(fsf0 + 99);

    const auto *fsd_51 = buffer.data(fsd + 51);
    const auto *fsd_57 = buffer.data(fsd + 57);
    const auto *fsd_59 = buffer.data(fsd + 59);

    const auto *fsf1_96 = buffer.data(fsf1 + 96);
    const auto *fsf1_99 = buffer.data(fsf1 + 99);

    const auto *gsp0_42 = buffer.data(gsp0 + 42);
    const auto *gsp0_43 = buffer.data(gsp0 + 43);
    const auto *gsp0_44 = buffer.data(gsp0 + 44);

    const auto *gsp1_42 = buffer.data(gsp1 + 42);
    const auto *gsp1_43 = buffer.data(gsp1 + 43);
    const auto *gsp1_44 = buffer.data(gsp1 + 44);

    const auto *gsd_81 = buffer.data(gsd + 81);
    const auto *gsd_83 = buffer.data(gsd + 83);
    const auto *gsd_84 = buffer.data(gsd + 84);
    const auto *gsd_86 = buffer.data(gsd + 86);
    const auto *gsd_87 = buffer.data(gsd + 87);
    const auto *gsd_88 = buffer.data(gsd + 88);
    const auto *gsd_89 = buffer.data(gsd + 89);

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pa_y, pc_y, pc_z, fsf0_96, fsf0_99, \
                         fsd_51, fsd_57, fsd_59, fsf1_96, fsf1_99, gsd_81, \
                         gsd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = pa_y[k] * fsf0_96[k]
                   + f_6 * fsd_57[k]
                   - f_4 * pc_y[k] * fsf1_96[k];

        t_137[k] = f_6 * fsd_51[k]
                   + f_3 * pc_z[k] * gsd_81[k];

        t_138[k] = f_5 * fsd_59[k]
                   + f_3 * pc_y[k] * gsd_83[k];

        t_139[k] = pa_y[k] * fsf0_99[k]
                   - f_4 * pc_y[k] * fsf1_99[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, pc_x, pc_y, gsp0_42, gsp0_44, \
                         gsp1_42, gsp1_44, gsd_84, gsd_86, gsd_87, \
                         gsd_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_1 * gsp0_42[k]
                   - f_2 * gsp1_42[k]
                   + f_3 * pc_x[k] * gsd_84[k];

        t_141[k] = f_3 * pc_y[k] * gsd_84[k];

        t_142[k] = f_7 * gsp0_44[k]
                   - f_8 * gsp1_44[k]
                   + f_3 * pc_x[k] * gsd_86[k];

        t_143[k] = f_3 * pc_x[k] * gsd_87[k];

        t_144[k] = f_3 * pc_x[k] * gsd_88[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, pc_x, pc_y, pc_z, fsd_59, gsp0_43, \
                         gsp0_44, gsp1_43, gsp1_44, gsd_87, gsd_88, \
                         gsd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_3 * pc_x[k] * gsd_89[k];

        t_146[k] = f_1 * gsp0_43[k]
                   - f_2 * gsp1_43[k]
                   + f_3 * pc_y[k] * gsd_87[k];

        t_147[k] = f_7 * gsp0_44[k]
                   - f_8 * gsp1_44[k]
                   + f_3 * pc_y[k] * gsd_88[k];

        t_148[k] = f_3 * pc_y[k] * gsd_89[k];

        t_149[k] = f_0 * fsd_59[k]
                   + f_1 * gsp0_44[k]
                   - f_2 * gsp1_44[k]
                   + f_3 * pc_z[k] * gsd_89[k];
    }
}

auto
compute_prim_gsf_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t fsf0, const size_t fsd,
                                                   const size_t fsf1, const size_t gsp0,
                                                   const size_t gsp1, const size_t gsd,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_gsf_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, fsf0, fsd,
                                                              fsf1, gsp0, gsp1, gsd, ncols,
                                                              gamma, p, q);

    compute_prim_gsf_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, fsf0, fsd,
                                                              fsf1, gsp0, gsp1, gsd, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
