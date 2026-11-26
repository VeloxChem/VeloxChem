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


#include "SimdThreeCenterElectronRepulsionVrrRecFSG.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_fsg_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t dsg0,
                                                          const size_t dsf, const size_t dsg1,
                                                          const size_t fsd0, const size_t fsd1,
                                                          const size_t fsf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 1.0 / gamma;
    const auto f_10 = p / (gamma * q);
    const auto f_11 = 2.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dsg0_0 = buffer.data(dsg0 + 0);
    const auto *dsg0_3 = buffer.data(dsg0 + 3);
    const auto *dsg0_5 = buffer.data(dsg0 + 5);
    const auto *dsg0_10 = buffer.data(dsg0 + 10);
    const auto *dsg0_14 = buffer.data(dsg0 + 14);
    const auto *dsg0_18 = buffer.data(dsg0 + 18);
    const auto *dsg0_30 = buffer.data(dsg0 + 30);
    const auto *dsg0_35 = buffer.data(dsg0 + 35);
    const auto *dsg0_45 = buffer.data(dsg0 + 45);
    const auto *dsg0_46 = buffer.data(dsg0 + 46);
    const auto *dsg0_48 = buffer.data(dsg0 + 48);
    const auto *dsg0_55 = buffer.data(dsg0 + 55);
    const auto *dsg0_57 = buffer.data(dsg0 + 57);
    const auto *dsg0_59 = buffer.data(dsg0 + 59);
    const auto *dsg0_70 = buffer.data(dsg0 + 70);
    const auto *dsg0_72 = buffer.data(dsg0 + 72);
    const auto *dsg0_74 = buffer.data(dsg0 + 74);
    const auto *dsg0_75 = buffer.data(dsg0 + 75);
    const auto *dsg0_77 = buffer.data(dsg0 + 77);
    const auto *dsg0_80 = buffer.data(dsg0 + 80);
    const auto *dsg0_85 = buffer.data(dsg0 + 85);
    const auto *dsg0_86 = buffer.data(dsg0 + 86);
    const auto *dsg0_87 = buffer.data(dsg0 + 87);
    const auto *dsg0_89 = buffer.data(dsg0 + 89);

    const auto *dsf_0 = buffer.data(dsf + 0);
    const auto *dsf_1 = buffer.data(dsf + 1);
    const auto *dsf_2 = buffer.data(dsf + 2);
    const auto *dsf_6 = buffer.data(dsf + 6);
    const auto *dsf_9 = buffer.data(dsf + 9);
    const auto *dsf_10 = buffer.data(dsf + 10);
    const auto *dsf_16 = buffer.data(dsf + 16);
    const auto *dsf_18 = buffer.data(dsf + 18);
    const auto *dsf_19 = buffer.data(dsf + 19);
    const auto *dsf_20 = buffer.data(dsf + 20);
    const auto *dsf_22 = buffer.data(dsf + 22);
    const auto *dsf_26 = buffer.data(dsf + 26);
    const auto *dsf_27 = buffer.data(dsf + 27);
    const auto *dsf_29 = buffer.data(dsf + 29);
    const auto *dsf_30 = buffer.data(dsf + 30);
    const auto *dsf_33 = buffer.data(dsf + 33);
    const auto *dsf_36 = buffer.data(dsf + 36);
    const auto *dsf_37 = buffer.data(dsf + 37);
    const auto *dsf_38 = buffer.data(dsf + 38);
    const auto *dsf_39 = buffer.data(dsf + 39);
    const auto *dsf_46 = buffer.data(dsf + 46);
    const auto *dsf_47 = buffer.data(dsf + 47);
    const auto *dsf_48 = buffer.data(dsf + 48);
    const auto *dsf_49 = buffer.data(dsf + 49);
    const auto *dsf_50 = buffer.data(dsf + 50);
    const auto *dsf_55 = buffer.data(dsf + 55);
    const auto *dsf_56 = buffer.data(dsf + 56);
    const auto *dsf_57 = buffer.data(dsf + 57);
    const auto *dsf_58 = buffer.data(dsf + 58);
    const auto *dsf_59 = buffer.data(dsf + 59);

    const auto *dsg1_0 = buffer.data(dsg1 + 0);
    const auto *dsg1_3 = buffer.data(dsg1 + 3);
    const auto *dsg1_5 = buffer.data(dsg1 + 5);
    const auto *dsg1_10 = buffer.data(dsg1 + 10);
    const auto *dsg1_14 = buffer.data(dsg1 + 14);
    const auto *dsg1_18 = buffer.data(dsg1 + 18);
    const auto *dsg1_30 = buffer.data(dsg1 + 30);
    const auto *dsg1_35 = buffer.data(dsg1 + 35);
    const auto *dsg1_45 = buffer.data(dsg1 + 45);
    const auto *dsg1_46 = buffer.data(dsg1 + 46);
    const auto *dsg1_48 = buffer.data(dsg1 + 48);
    const auto *dsg1_55 = buffer.data(dsg1 + 55);
    const auto *dsg1_57 = buffer.data(dsg1 + 57);
    const auto *dsg1_59 = buffer.data(dsg1 + 59);
    const auto *dsg1_70 = buffer.data(dsg1 + 70);
    const auto *dsg1_72 = buffer.data(dsg1 + 72);
    const auto *dsg1_74 = buffer.data(dsg1 + 74);
    const auto *dsg1_75 = buffer.data(dsg1 + 75);
    const auto *dsg1_77 = buffer.data(dsg1 + 77);
    const auto *dsg1_80 = buffer.data(dsg1 + 80);
    const auto *dsg1_85 = buffer.data(dsg1 + 85);
    const auto *dsg1_86 = buffer.data(dsg1 + 86);
    const auto *dsg1_87 = buffer.data(dsg1 + 87);
    const auto *dsg1_89 = buffer.data(dsg1 + 89);

    const auto *fsd0_0 = buffer.data(fsd0 + 0);
    const auto *fsd0_3 = buffer.data(fsd0 + 3);
    const auto *fsd0_5 = buffer.data(fsd0 + 5);
    const auto *fsd0_9 = buffer.data(fsd0 + 9);
    const auto *fsd0_16 = buffer.data(fsd0 + 16);
    const auto *fsd0_17 = buffer.data(fsd0 + 17);
    const auto *fsd0_18 = buffer.data(fsd0 + 18);
    const auto *fsd0_30 = buffer.data(fsd0 + 30);
    const auto *fsd0_36 = buffer.data(fsd0 + 36);
    const auto *fsd0_37 = buffer.data(fsd0 + 37);
    const auto *fsd0_39 = buffer.data(fsd0 + 39);
    const auto *fsd0_41 = buffer.data(fsd0 + 41);
    const auto *fsd0_44 = buffer.data(fsd0 + 44);
    const auto *fsd0_46 = buffer.data(fsd0 + 46);
    const auto *fsd0_47 = buffer.data(fsd0 + 47);
    const auto *fsd0_49 = buffer.data(fsd0 + 49);
    const auto *fsd0_51 = buffer.data(fsd0 + 51);
    const auto *fsd0_52 = buffer.data(fsd0 + 52);

    const auto *fsd1_0 = buffer.data(fsd1 + 0);
    const auto *fsd1_3 = buffer.data(fsd1 + 3);
    const auto *fsd1_5 = buffer.data(fsd1 + 5);
    const auto *fsd1_9 = buffer.data(fsd1 + 9);
    const auto *fsd1_16 = buffer.data(fsd1 + 16);
    const auto *fsd1_17 = buffer.data(fsd1 + 17);
    const auto *fsd1_18 = buffer.data(fsd1 + 18);
    const auto *fsd1_30 = buffer.data(fsd1 + 30);
    const auto *fsd1_36 = buffer.data(fsd1 + 36);
    const auto *fsd1_37 = buffer.data(fsd1 + 37);
    const auto *fsd1_39 = buffer.data(fsd1 + 39);
    const auto *fsd1_41 = buffer.data(fsd1 + 41);
    const auto *fsd1_44 = buffer.data(fsd1 + 44);
    const auto *fsd1_46 = buffer.data(fsd1 + 46);
    const auto *fsd1_47 = buffer.data(fsd1 + 47);
    const auto *fsd1_49 = buffer.data(fsd1 + 49);
    const auto *fsd1_51 = buffer.data(fsd1 + 51);
    const auto *fsd1_52 = buffer.data(fsd1 + 52);

    const auto *fsf_0 = buffer.data(fsf + 0);
    const auto *fsf_1 = buffer.data(fsf + 1);
    const auto *fsf_2 = buffer.data(fsf + 2);
    const auto *fsf_3 = buffer.data(fsf + 3);
    const auto *fsf_5 = buffer.data(fsf + 5);
    const auto *fsf_6 = buffer.data(fsf + 6);
    const auto *fsf_8 = buffer.data(fsf + 8);
    const auto *fsf_9 = buffer.data(fsf + 9);
    const auto *fsf_10 = buffer.data(fsf + 10);
    const auto *fsf_11 = buffer.data(fsf + 11);
    const auto *fsf_13 = buffer.data(fsf + 13);
    const auto *fsf_16 = buffer.data(fsf + 16);
    const auto *fsf_17 = buffer.data(fsf + 17);
    const auto *fsf_18 = buffer.data(fsf + 18);
    const auto *fsf_19 = buffer.data(fsf + 19);
    const auto *fsf_20 = buffer.data(fsf + 20);
    const auto *fsf_22 = buffer.data(fsf + 22);
    const auto *fsf_25 = buffer.data(fsf + 25);
    const auto *fsf_26 = buffer.data(fsf + 26);
    const auto *fsf_27 = buffer.data(fsf + 27);
    const auto *fsf_28 = buffer.data(fsf + 28);
    const auto *fsf_29 = buffer.data(fsf + 29);
    const auto *fsf_30 = buffer.data(fsf + 30);
    const auto *fsf_31 = buffer.data(fsf + 31);
    const auto *fsf_32 = buffer.data(fsf + 32);
    const auto *fsf_33 = buffer.data(fsf + 33);
    const auto *fsf_36 = buffer.data(fsf + 36);
    const auto *fsf_38 = buffer.data(fsf + 38);
    const auto *fsf_39 = buffer.data(fsf + 39);
    const auto *fsf_40 = buffer.data(fsf + 40);
    const auto *fsf_42 = buffer.data(fsf + 42);
    const auto *fsf_46 = buffer.data(fsf + 46);
    const auto *fsf_47 = buffer.data(fsf + 47);
    const auto *fsf_48 = buffer.data(fsf + 48);
    const auto *fsf_49 = buffer.data(fsf + 49);
    const auto *fsf_50 = buffer.data(fsf + 50);
    const auto *fsf_51 = buffer.data(fsf + 51);
    const auto *fsf_52 = buffer.data(fsf + 52);
    const auto *fsf_55 = buffer.data(fsf + 55);
    const auto *fsf_56 = buffer.data(fsf + 56);
    const auto *fsf_57 = buffer.data(fsf + 57);
    const auto *fsf_59 = buffer.data(fsf + 59);
    const auto *fsf_60 = buffer.data(fsf + 60);
    const auto *fsf_61 = buffer.data(fsf + 61);
    const auto *fsf_63 = buffer.data(fsf + 63);
    const auto *fsf_65 = buffer.data(fsf + 65);
    const auto *fsf_66 = buffer.data(fsf + 66);
    const auto *fsf_67 = buffer.data(fsf + 67);
    const auto *fsf_68 = buffer.data(fsf + 68);
    const auto *fsf_69 = buffer.data(fsf + 69);
    const auto *fsf_72 = buffer.data(fsf + 72);
    const auto *fsf_74 = buffer.data(fsf + 74);
    const auto *fsf_75 = buffer.data(fsf + 75);
    const auto *fsf_76 = buffer.data(fsf + 76);
    const auto *fsf_77 = buffer.data(fsf + 77);
    const auto *fsf_78 = buffer.data(fsf + 78);
    const auto *fsf_79 = buffer.data(fsf + 79);
    const auto *fsf_81 = buffer.data(fsf + 81);
    const auto *fsf_83 = buffer.data(fsf + 83);
    const auto *fsf_84 = buffer.data(fsf + 84);
    const auto *fsf_86 = buffer.data(fsf + 86);
    const auto *fsf_87 = buffer.data(fsf + 87);
    const auto *fsf_88 = buffer.data(fsf + 88);
    const auto *fsf_89 = buffer.data(fsf + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, dsf_0, fsd0_0, \
                         fsd1_0, fsf_0, fsf_1, fsf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dsf_0[k]
                 + f_1 * fsd0_0[k]
                 - f_2 * fsd1_0[k]
                 + f_3 * pc_x[k] * fsf_0[k];

        t_1[k] = f_3 * pc_y[k] * fsf_0[k];

        t_2[k] = f_3 * pc_z[k] * fsf_0[k];

        t_3[k] = f_4 * fsd0_0[k]
                 - f_5 * fsd1_0[k]
                 + f_3 * pc_y[k] * fsf_1[k];

        t_4[k] = f_3 * pc_y[k] * fsf_2[k];

        t_5[k] = f_4 * fsd0_0[k]
                 - f_5 * fsd1_0[k]
                 + f_3 * pc_z[k] * fsf_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_x, pc_y, pc_z, dsf_6, dsf_9, fsd0_3, \
                         fsd1_3, fsf_3, fsf_5, fsf_6, fsf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * dsf_6[k]
                 + f_3 * pc_x[k] * fsf_6[k];

        t_7[k] = f_3 * pc_z[k] * fsf_3[k];

        t_8[k] = f_3 * pc_y[k] * fsf_5[k];

        t_9[k] = f_0 * dsf_9[k]
                 + f_3 * pc_x[k] * fsf_9[k];

        t_10[k] = f_1 * fsd0_3[k]
                  - f_2 * fsd1_3[k]
                  + f_3 * pc_y[k] * fsf_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pc_y, pc_z, dsg0_0, dsg1_0, \
                         fsd0_5, fsd1_5, fsf_6, fsf_8, fsf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * fsf_6[k];

        t_12[k] = f_4 * fsd0_5[k]
                  - f_5 * fsd1_5[k]
                  + f_3 * pc_y[k] * fsf_8[k];

        t_13[k] = f_3 * pc_y[k] * fsf_9[k];

        t_14[k] = f_1 * fsd0_5[k]
                  - f_2 * fsd1_5[k]
                  + f_3 * pc_z[k] * fsf_9[k];

        t_15[k] = pa_y[k] * dsg0_0[k]
                  - f_6 * pc_y[k] * dsg1_0[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_y, pc_y, pc_z, dsg0_3, dsg0_5, \
                         dsf_0, dsf_1, dsg1_3, dsg1_5, fsf_10, fsf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_7 * dsf_0[k]
                  + f_3 * pc_y[k] * fsf_10[k];

        t_17[k] = f_3 * pc_z[k] * fsf_10[k];

        t_18[k] = pa_y[k] * dsg0_3[k]
                  + f_8 * dsf_1[k]
                  - f_6 * pc_y[k] * dsg1_3[k];

        t_19[k] = f_3 * pc_z[k] * fsf_11[k];

        t_20[k] = pa_y[k] * dsg0_5[k]
                  - f_6 * pc_y[k] * dsg1_5[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_x, pc_z, dsf_16, dsf_18, dsf_19, fsf_13, \
                         fsf_16, fsf_18, fsf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_8 * dsf_16[k]
                  + f_3 * pc_x[k] * fsf_16[k];

        t_22[k] = f_3 * pc_z[k] * fsf_13[k];

        t_23[k] = f_8 * dsf_18[k]
                  + f_3 * pc_x[k] * fsf_18[k];

        t_24[k] = f_8 * dsf_19[k]
                  + f_3 * pc_x[k] * fsf_19[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pc_y, pc_z, dsf_6, dsf_9, fsd0_9, fsd1_9, \
                         fsf_16, fsf_17, fsf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_7 * dsf_6[k]
                  + f_1 * fsd0_9[k]
                  - f_2 * fsd1_9[k]
                  + f_3 * pc_y[k] * fsf_16[k];

        t_26[k] = f_3 * pc_z[k] * fsf_16[k];

        t_27[k] = f_4 * fsd0_9[k]
                  - f_5 * fsd1_9[k]
                  + f_3 * pc_z[k] * fsf_17[k];

        t_28[k] = f_7 * dsf_9[k]
                  + f_3 * pc_y[k] * fsf_19[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_y, pa_z, pc_y, pc_z, dsg0_0, dsg0_14, \
                         dsf_0, dsg1_0, dsg1_14, fsf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_y[k] * dsg0_14[k]
                  - f_6 * pc_y[k] * dsg1_14[k];

        t_30[k] = pa_z[k] * dsg0_0[k]
                  - f_6 * pc_z[k] * dsg1_0[k];

        t_31[k] = f_3 * pc_y[k] * fsf_20[k];

        t_32[k] = f_7 * dsf_0[k]
                  + f_3 * pc_z[k] * fsf_20[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_z, pc_x, pc_y, pc_z, dsg0_3, dsg0_5, \
                         dsf_2, dsf_26, dsg1_3, dsg1_5, fsf_22, \
                         fsf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_z[k] * dsg0_3[k]
                  - f_6 * pc_z[k] * dsg1_3[k];

        t_34[k] = f_3 * pc_y[k] * fsf_22[k];

        t_35[k] = pa_z[k] * dsg0_5[k]
                  + f_8 * dsf_2[k]
                  - f_6 * pc_z[k] * dsg1_5[k];

        t_36[k] = f_8 * dsf_26[k]
                  + f_3 * pc_x[k] * fsf_26[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_z, pc_x, pc_y, pc_z, dsg0_10, dsf_27, \
                         dsf_29, dsg1_10, fsf_25, fsf_27, fsf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_8 * dsf_27[k]
                  + f_3 * pc_x[k] * fsf_27[k];

        t_38[k] = f_3 * pc_y[k] * fsf_25[k];

        t_39[k] = f_8 * dsf_29[k]
                  + f_3 * pc_x[k] * fsf_29[k];

        t_40[k] = pa_z[k] * dsg0_10[k]
                  - f_6 * pc_z[k] * dsg1_10[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pc_y, pc_z, dsf_9, fsd0_16, fsd0_17, fsd1_16, \
                         fsd1_17, fsf_27, fsf_28, fsf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_9 * fsd0_16[k]
                  - f_10 * fsd1_16[k]
                  + f_3 * pc_y[k] * fsf_27[k];

        t_42[k] = f_4 * fsd0_17[k]
                  - f_5 * fsd1_17[k]
                  + f_3 * pc_y[k] * fsf_28[k];

        t_43[k] = f_3 * pc_y[k] * fsf_29[k];

        t_44[k] = f_7 * dsf_9[k]
                  + f_1 * fsd0_17[k]
                  - f_2 * fsd1_17[k]
                  + f_3 * pc_z[k] * fsf_29[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pa_x, pc_x, pc_y, pc_z, dsg0_45, dsg0_48, \
                         dsf_10, dsf_30, dsf_33, dsg1_45, dsg1_48, \
                         fsf_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pa_x[k] * dsg0_45[k]
                  + f_11 * dsf_30[k]
                  - f_6 * pc_x[k] * dsg1_45[k];

        t_46[k] = f_8 * dsf_10[k]
                  + f_3 * pc_y[k] * fsf_30[k];

        t_47[k] = f_3 * pc_z[k] * fsf_30[k];

        t_48[k] = pa_x[k] * dsg0_48[k]
                  + f_8 * dsf_33[k]
                  - f_6 * pc_x[k] * dsg1_48[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, pc_x, pc_z, dsf_36, dsf_38, fsd0_18, \
                         fsd1_18, fsf_31, fsf_32, fsf_33, fsf_36, \
                         fsf_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_3 * pc_z[k] * fsf_31[k];

        t_50[k] = f_4 * fsd0_18[k]
                  - f_5 * fsd1_18[k]
                  + f_3 * pc_z[k] * fsf_32[k];

        t_51[k] = f_7 * dsf_36[k]
                  + f_3 * pc_x[k] * fsf_36[k];

        t_52[k] = f_3 * pc_z[k] * fsf_33[k];

        t_53[k] = f_7 * dsf_38[k]
                  + f_3 * pc_x[k] * fsf_38[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_x, pc_x, pc_z, dsg0_55, dsg0_57, dsf_39, \
                         dsg1_55, dsg1_57, fsf_36, fsf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_7 * dsf_39[k]
                  + f_3 * pc_x[k] * fsf_39[k];

        t_55[k] = pa_x[k] * dsg0_55[k]
                  - f_6 * pc_x[k] * dsg1_55[k];

        t_56[k] = f_3 * pc_z[k] * fsf_36[k];

        t_57[k] = pa_x[k] * dsg0_57[k]
                  - f_6 * pc_x[k] * dsg1_57[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pa_x, pa_y, pc_x, pc_y, dsg0_30, dsg0_59, \
                         dsf_19, dsf_20, dsg1_30, dsg1_59, fsf_39, \
                         fsf_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_8 * dsf_19[k]
                  + f_3 * pc_y[k] * fsf_39[k];

        t_59[k] = pa_x[k] * dsg0_59[k]
                  - f_6 * pc_x[k] * dsg1_59[k];

        t_60[k] = pa_y[k] * dsg0_30[k]
                  - f_6 * pc_y[k] * dsg1_30[k];

        t_61[k] = f_7 * dsf_20[k]
                  + f_3 * pc_y[k] * fsf_40[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_y, pa_z, pc_y, pc_z, dsg0_18, dsg0_35, \
                         dsf_10, dsf_22, dsg1_18, dsg1_35, fsf_40, \
                         fsf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_7 * dsf_10[k]
                  + f_3 * pc_z[k] * fsf_40[k];

        t_63[k] = pa_z[k] * dsg0_18[k]
                  - f_6 * pc_z[k] * dsg1_18[k];

        t_64[k] = f_7 * dsf_22[k]
                  + f_3 * pc_y[k] * fsf_42[k];

        t_65[k] = pa_y[k] * dsg0_35[k]
                  - f_6 * pc_y[k] * dsg1_35[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pc_x, dsf_46, dsf_47, dsf_48, dsf_49, fsf_46, \
                         fsf_47, fsf_48, fsf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_7 * dsf_46[k]
                  + f_3 * pc_x[k] * fsf_46[k];

        t_67[k] = f_7 * dsf_47[k]
                  + f_3 * pc_x[k] * fsf_47[k];

        t_68[k] = f_7 * dsf_48[k]
                  + f_3 * pc_x[k] * fsf_48[k];

        t_69[k] = f_7 * dsf_49[k]
                  + f_3 * pc_x[k] * fsf_49[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pa_x, pc_x, pc_y, pc_z, dsg0_70, dsg0_72, \
                         dsf_16, dsf_29, dsg1_70, dsg1_72, fsf_46, \
                         fsf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pa_x[k] * dsg0_70[k]
                  - f_6 * pc_x[k] * dsg1_70[k];

        t_71[k] = f_7 * dsf_16[k]
                  + f_3 * pc_z[k] * fsf_46[k];

        t_72[k] = pa_x[k] * dsg0_72[k]
                  - f_6 * pc_x[k] * dsg1_72[k];

        t_73[k] = f_7 * dsf_29[k]
                  + f_3 * pc_y[k] * fsf_49[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pa_x, pc_x, pc_y, pc_z, dsg0_74, dsg0_75, \
                         dsf_20, dsf_50, dsg1_74, dsg1_75, fsf_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = pa_x[k] * dsg0_74[k]
                  - f_6 * pc_x[k] * dsg1_74[k];

        t_75[k] = pa_x[k] * dsg0_75[k]
                  + f_11 * dsf_50[k]
                  - f_6 * pc_x[k] * dsg1_75[k];

        t_76[k] = f_3 * pc_y[k] * fsf_50[k];

        t_77[k] = f_8 * dsf_20[k]
                  + f_3 * pc_z[k] * fsf_50[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_x, pc_x, pc_y, dsg0_80, dsf_55, dsf_56, \
                         dsg1_80, fsd0_30, fsd1_30, fsf_51, fsf_52, \
                         fsf_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_4 * fsd0_30[k]
                  - f_5 * fsd1_30[k]
                  + f_3 * pc_y[k] * fsf_51[k];

        t_79[k] = f_3 * pc_y[k] * fsf_52[k];

        t_80[k] = pa_x[k] * dsg0_80[k]
                  + f_8 * dsf_55[k]
                  - f_6 * pc_x[k] * dsg1_80[k];

        t_81[k] = f_7 * dsf_56[k]
                  + f_3 * pc_x[k] * fsf_56[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pa_x, pc_x, pc_y, dsg0_85, dsf_57, dsf_59, \
                         dsg1_85, fsf_55, fsf_57, fsf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_7 * dsf_57[k]
                  + f_3 * pc_x[k] * fsf_57[k];

        t_83[k] = f_3 * pc_y[k] * fsf_55[k];

        t_84[k] = f_7 * dsf_59[k]
                  + f_3 * pc_x[k] * fsf_59[k];

        t_85[k] = pa_x[k] * dsg0_85[k]
                  - f_6 * pc_x[k] * dsg1_85[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pa_x, pc_x, pc_y, dsg0_86, dsg0_87, dsg0_89, \
                         dsg1_86, dsg1_87, dsg1_89, fsf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pa_x[k] * dsg0_86[k]
                  - f_6 * pc_x[k] * dsg1_86[k];

        t_87[k] = pa_x[k] * dsg0_87[k]
                  - f_6 * pc_x[k] * dsg1_87[k];

        t_88[k] = f_3 * pc_y[k] * fsf_59[k];

        t_89[k] = pa_x[k] * dsg0_89[k]
                  - f_6 * pc_x[k] * dsg1_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, pc_x, pc_z, fsd0_36, fsd0_37, fsd0_39, \
                         fsd1_36, fsd1_37, fsd1_39, fsf_60, fsf_61, \
                         fsf_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_1 * fsd0_36[k]
                  - f_2 * fsd1_36[k]
                  + f_3 * pc_x[k] * fsf_60[k];

        t_91[k] = f_9 * fsd0_37[k]
                  - f_10 * fsd1_37[k]
                  + f_3 * pc_x[k] * fsf_61[k];

        t_92[k] = f_3 * pc_z[k] * fsf_60[k];

        t_93[k] = f_4 * fsd0_39[k]
                  - f_5 * fsd1_39[k]
                  + f_3 * pc_x[k] * fsf_63[k];

        t_94[k] = f_3 * pc_z[k] * fsf_61[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, pc_x, fsd0_41, fsd1_41, fsf_65, fsf_66, \
                         fsf_67, fsf_68, fsf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_4 * fsd0_41[k]
                  - f_5 * fsd1_41[k]
                  + f_3 * pc_x[k] * fsf_65[k];

        t_96[k] = f_3 * pc_x[k] * fsf_66[k];

        t_97[k] = f_3 * pc_x[k] * fsf_67[k];

        t_98[k] = f_3 * pc_x[k] * fsf_68[k];

        t_99[k] = f_3 * pc_x[k] * fsf_69[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, t_104, pc_y, pc_z, dsf_36, dsf_39, \
                         fsd0_39, fsd0_41, fsd1_39, fsd1_41, fsf_66, fsf_67, \
                         fsf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_0 * dsf_36[k]
                   + f_1 * fsd0_39[k]
                   - f_2 * fsd1_39[k]
                   + f_3 * pc_y[k] * fsf_66[k];

        t_101[k] = f_3 * pc_z[k] * fsf_66[k];

        t_102[k] = f_4 * fsd0_39[k]
                   - f_5 * fsd1_39[k]
                   + f_3 * pc_z[k] * fsf_67[k];

        t_103[k] = f_0 * dsf_39[k]
                   + f_3 * pc_y[k] * fsf_69[k];

        t_104[k] = f_1 * fsd0_41[k]
                   - f_2 * fsd1_41[k]
                   + f_3 * pc_z[k] * fsf_69[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pa_z, pc_x, pc_z, dsg0_45, dsg0_46, \
                         dsg0_48, dsg1_45, dsg1_46, dsg1_48, fsd0_44, fsd1_44, \
                         fsf_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = pa_z[k] * dsg0_45[k]
                   - f_6 * pc_z[k] * dsg1_45[k];

        t_106[k] = pa_z[k] * dsg0_46[k]
                   - f_6 * pc_z[k] * dsg1_46[k];

        t_107[k] = f_9 * fsd0_44[k]
                   - f_10 * fsd1_44[k]
                   + f_3 * pc_x[k] * fsf_72[k];

        t_108[k] = pa_z[k] * dsg0_48[k]
                   - f_6 * pc_z[k] * dsg1_48[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, t_113, pc_x, fsd0_46, fsd0_47, fsd1_46, \
                         fsd1_47, fsf_74, fsf_75, fsf_76, fsf_77, \
                         fsf_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_4 * fsd0_46[k]
                   - f_5 * fsd1_46[k]
                   + f_3 * pc_x[k] * fsf_74[k];

        t_110[k] = f_4 * fsd0_47[k]
                   - f_5 * fsd1_47[k]
                   + f_3 * pc_x[k] * fsf_75[k];

        t_111[k] = f_3 * pc_x[k] * fsf_76[k];

        t_112[k] = f_3 * pc_x[k] * fsf_77[k];

        t_113[k] = f_3 * pc_x[k] * fsf_78[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pa_z, pc_x, pc_z, dsg0_55, dsg0_57, \
                         dsf_36, dsf_37, dsg1_55, dsg1_57, fsf_76, \
                         fsf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_3 * pc_x[k] * fsf_79[k];

        t_115[k] = pa_z[k] * dsg0_55[k]
                   - f_6 * pc_z[k] * dsg1_55[k];

        t_116[k] = f_7 * dsf_36[k]
                   + f_3 * pc_z[k] * fsf_76[k];

        t_117[k] = pa_z[k] * dsg0_57[k]
                   + f_8 * dsf_37[k]
                   - f_6 * pc_z[k] * dsg1_57[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pa_y, pc_y, pc_z, dsg0_75, dsf_39, dsf_49, \
                         dsg1_75, fsd0_47, fsd1_47, fsf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_8 * dsf_49[k]
                   + f_3 * pc_y[k] * fsf_79[k];

        t_119[k] = f_7 * dsf_39[k]
                   + f_1 * fsd0_47[k]
                   - f_2 * fsd1_47[k]
                   + f_3 * pc_z[k] * fsf_79[k];

        t_120[k] = pa_y[k] * dsg0_75[k]
                   - f_6 * pc_y[k] * dsg1_75[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, pa_y, pc_x, pc_y, dsg0_77, dsg1_77, fsd0_49, \
                         fsd0_51, fsd1_49, fsd1_51, fsf_81, fsf_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_9 * fsd0_49[k]
                   - f_10 * fsd1_49[k]
                   + f_3 * pc_x[k] * fsf_81[k];

        t_122[k] = pa_y[k] * dsg0_77[k]
                   - f_6 * pc_y[k] * dsg1_77[k];

        t_123[k] = f_4 * fsd0_51[k]
                   - f_5 * fsd1_51[k]
                   + f_3 * pc_x[k] * fsf_83[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, pa_y, pc_x, pc_y, dsg0_80, \
                         dsg1_80, fsd0_52, fsd1_52, fsf_84, fsf_86, fsf_87, \
                         fsf_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_4 * fsd0_52[k]
                   - f_5 * fsd1_52[k]
                   + f_3 * pc_x[k] * fsf_84[k];

        t_125[k] = pa_y[k] * dsg0_80[k]
                   - f_6 * pc_y[k] * dsg1_80[k];

        t_126[k] = f_3 * pc_x[k] * fsf_86[k];

        t_127[k] = f_3 * pc_x[k] * fsf_87[k];

        t_128[k] = f_3 * pc_x[k] * fsf_88[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, pa_y, pc_x, pc_y, pc_z, dsg0_85, dsf_46, dsf_56, \
                         dsg1_85, fsf_86, fsf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_3 * pc_x[k] * fsf_89[k];

        t_130[k] = pa_y[k] * dsg0_85[k]
                   + f_11 * dsf_56[k]
                   - f_6 * pc_y[k] * dsg1_85[k];

        t_131[k] = f_8 * dsf_46[k]
                   + f_3 * pc_z[k] * fsf_86[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, pa_y, pc_y, dsg0_87, dsg0_89, dsf_58, dsf_59, \
                         dsg1_87, dsg1_89, fsf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = pa_y[k] * dsg0_87[k]
                   + f_8 * dsf_58[k]
                   - f_6 * pc_y[k] * dsg1_87[k];

        t_133[k] = f_7 * dsf_59[k]
                   + f_3 * pc_y[k] * fsf_89[k];

        t_134[k] = pa_y[k] * dsg0_89[k]
                   - f_6 * pc_y[k] * dsg1_89[k];
    }
}

static auto
compute_prim_fsg_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pc,
                                                          const size_t dsf, const size_t fsd0,
                                                          const size_t fsd1, const size_t fsf,
                                                          const size_t ncols, const double gamma,
                                                          const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_9 = 1.0 / gamma;
    const auto f_10 = p / (gamma * q);

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dsf_59 = buffer.data(dsf + 59);

    const auto *fsd0_54 = buffer.data(fsd0 + 54);
    const auto *fsd0_56 = buffer.data(fsd0 + 56);
    const auto *fsd0_57 = buffer.data(fsd0 + 57);
    const auto *fsd0_58 = buffer.data(fsd0 + 58);
    const auto *fsd0_59 = buffer.data(fsd0 + 59);

    const auto *fsd1_54 = buffer.data(fsd1 + 54);
    const auto *fsd1_56 = buffer.data(fsd1 + 56);
    const auto *fsd1_57 = buffer.data(fsd1 + 57);
    const auto *fsd1_58 = buffer.data(fsd1 + 58);
    const auto *fsd1_59 = buffer.data(fsd1 + 59);

    const auto *fsf_90 = buffer.data(fsf + 90);
    const auto *fsf_92 = buffer.data(fsf + 92);
    const auto *fsf_93 = buffer.data(fsf + 93);
    const auto *fsf_95 = buffer.data(fsf + 95);
    const auto *fsf_96 = buffer.data(fsf + 96);
    const auto *fsf_97 = buffer.data(fsf + 97);
    const auto *fsf_98 = buffer.data(fsf + 98);
    const auto *fsf_99 = buffer.data(fsf + 99);

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, pc_x, pc_y, fsd0_54, fsd0_56, \
                         fsd0_57, fsd1_54, fsd1_56, fsd1_57, fsf_90, fsf_92, \
                         fsf_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_1 * fsd0_54[k]
                   - f_2 * fsd1_54[k]
                   + f_3 * pc_x[k] * fsf_90[k];

        t_136[k] = f_3 * pc_y[k] * fsf_90[k];

        t_137[k] = f_9 * fsd0_56[k]
                   - f_10 * fsd1_56[k]
                   + f_3 * pc_x[k] * fsf_92[k];

        t_138[k] = f_4 * fsd0_57[k]
                   - f_5 * fsd1_57[k]
                   + f_3 * pc_x[k] * fsf_93[k];

        t_139[k] = f_3 * pc_y[k] * fsf_92[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, pc_x, fsd0_59, fsd1_59, fsf_95, \
                         fsf_96, fsf_97, fsf_98, fsf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_4 * fsd0_59[k]
                   - f_5 * fsd1_59[k]
                   + f_3 * pc_x[k] * fsf_95[k];

        t_141[k] = f_3 * pc_x[k] * fsf_96[k];

        t_142[k] = f_3 * pc_x[k] * fsf_97[k];

        t_143[k] = f_3 * pc_x[k] * fsf_98[k];

        t_144[k] = f_3 * pc_x[k] * fsf_99[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, pc_y, fsd0_57, fsd0_58, fsd0_59, fsd1_57, \
                         fsd1_58, fsd1_59, fsf_96, fsf_97, fsf_98, \
                         fsf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_1 * fsd0_57[k]
                   - f_2 * fsd1_57[k]
                   + f_3 * pc_y[k] * fsf_96[k];

        t_146[k] = f_9 * fsd0_58[k]
                   - f_10 * fsd1_58[k]
                   + f_3 * pc_y[k] * fsf_97[k];

        t_147[k] = f_4 * fsd0_59[k]
                   - f_5 * fsd1_59[k]
                   + f_3 * pc_y[k] * fsf_98[k];

        t_148[k] = f_3 * pc_y[k] * fsf_99[k];
    }

#pragma omp simd aligned(t_149, pc_z, dsf_59, fsd0_59, fsd1_59, \
                         fsf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_0 * dsf_59[k]
                   + f_1 * fsd0_59[k]
                   - f_2 * fsd1_59[k]
                   + f_3 * pc_z[k] * fsf_99[k];
    }
}

auto
compute_prim_fsg_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t dsg0, const size_t dsf,
                                                   const size_t dsg1, const size_t fsd0,
                                                   const size_t fsd1, const size_t fsf,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_fsg_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, dsg0, dsf,
                                                              dsg1, fsd0, fsd1, fsf, ncols,
                                                              gamma, p, q);

    compute_prim_fsg_three_center_electron_repulsion_0_piece1(buffer, target, pc, dsf, fsd0,
                                                              fsd1, fsf, ncols, gamma, p, q);
}

}  // namespace simdt3ceri
