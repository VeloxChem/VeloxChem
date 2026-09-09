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


#include "SimdThreeCenterElectronRepulsionVrrRecGSG.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_gsg_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t fsg0,
                                                          const size_t fsf, const size_t fsg1,
                                                          const size_t gsd0, const size_t gsd1,
                                                          const size_t gsf, const size_t ncols,
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
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *fsg0_0 = buffer.data(fsg0 + 0);
    const auto *fsg0_3 = buffer.data(fsg0 + 3);
    const auto *fsg0_5 = buffer.data(fsg0 + 5);
    const auto *fsg0_10 = buffer.data(fsg0 + 10);
    const auto *fsg0_14 = buffer.data(fsg0 + 14);
    const auto *fsg0_18 = buffer.data(fsg0 + 18);
    const auto *fsg0_25 = buffer.data(fsg0 + 25);
    const auto *fsg0_30 = buffer.data(fsg0 + 30);
    const auto *fsg0_35 = buffer.data(fsg0 + 35);
    const auto *fsg0_44 = buffer.data(fsg0 + 44);
    const auto *fsg0_45 = buffer.data(fsg0 + 45);
    const auto *fsg0_48 = buffer.data(fsg0 + 48);
    const auto *fsg0_75 = buffer.data(fsg0 + 75);
    const auto *fsg0_80 = buffer.data(fsg0 + 80);
    const auto *fsg0_90 = buffer.data(fsg0 + 90);
    const auto *fsg0_93 = buffer.data(fsg0 + 93);
    const auto *fsg0_100 = buffer.data(fsg0 + 100);
    const auto *fsg0_102 = buffer.data(fsg0 + 102);
    const auto *fsg0_104 = buffer.data(fsg0 + 104);
    const auto *fsg0_110 = buffer.data(fsg0 + 110);
    const auto *fsg0_115 = buffer.data(fsg0 + 115);
    const auto *fsg0_117 = buffer.data(fsg0 + 117);
    const auto *fsg0_119 = buffer.data(fsg0 + 119);
    const auto *fsg0_123 = buffer.data(fsg0 + 123);
    const auto *fsg0_130 = buffer.data(fsg0 + 130);
    const auto *fsg0_132 = buffer.data(fsg0 + 132);

    const auto *fsf_0 = buffer.data(fsf + 0);
    const auto *fsf_1 = buffer.data(fsf + 1);
    const auto *fsf_2 = buffer.data(fsf + 2);
    const auto *fsf_6 = buffer.data(fsf + 6);
    const auto *fsf_9 = buffer.data(fsf + 9);
    const auto *fsf_10 = buffer.data(fsf + 10);
    const auto *fsf_16 = buffer.data(fsf + 16);
    const auto *fsf_18 = buffer.data(fsf + 18);
    const auto *fsf_19 = buffer.data(fsf + 19);
    const auto *fsf_20 = buffer.data(fsf + 20);
    const auto *fsf_22 = buffer.data(fsf + 22);
    const auto *fsf_26 = buffer.data(fsf + 26);
    const auto *fsf_27 = buffer.data(fsf + 27);
    const auto *fsf_28 = buffer.data(fsf + 28);
    const auto *fsf_29 = buffer.data(fsf + 29);
    const auto *fsf_30 = buffer.data(fsf + 30);
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
    const auto *fsf_52 = buffer.data(fsf + 52);
    const auto *fsf_55 = buffer.data(fsf + 55);
    const auto *fsf_56 = buffer.data(fsf + 56);
    const auto *fsf_57 = buffer.data(fsf + 57);
    const auto *fsf_59 = buffer.data(fsf + 59);
    const auto *fsf_60 = buffer.data(fsf + 60);
    const auto *fsf_63 = buffer.data(fsf + 63);
    const auto *fsf_66 = buffer.data(fsf + 66);
    const auto *fsf_68 = buffer.data(fsf + 68);
    const auto *fsf_69 = buffer.data(fsf + 69);
    const auto *fsf_75 = buffer.data(fsf + 75);
    const auto *fsf_76 = buffer.data(fsf + 76);
    const auto *fsf_77 = buffer.data(fsf + 77);
    const auto *fsf_78 = buffer.data(fsf + 78);
    const auto *fsf_79 = buffer.data(fsf + 79);
    const auto *fsf_83 = buffer.data(fsf + 83);
    const auto *fsf_86 = buffer.data(fsf + 86);
    const auto *fsf_87 = buffer.data(fsf + 87);
    const auto *fsf_88 = buffer.data(fsf + 88);
    const auto *fsf_89 = buffer.data(fsf + 89);

    const auto *fsg1_0 = buffer.data(fsg1 + 0);
    const auto *fsg1_3 = buffer.data(fsg1 + 3);
    const auto *fsg1_5 = buffer.data(fsg1 + 5);
    const auto *fsg1_10 = buffer.data(fsg1 + 10);
    const auto *fsg1_14 = buffer.data(fsg1 + 14);
    const auto *fsg1_18 = buffer.data(fsg1 + 18);
    const auto *fsg1_25 = buffer.data(fsg1 + 25);
    const auto *fsg1_30 = buffer.data(fsg1 + 30);
    const auto *fsg1_35 = buffer.data(fsg1 + 35);
    const auto *fsg1_44 = buffer.data(fsg1 + 44);
    const auto *fsg1_45 = buffer.data(fsg1 + 45);
    const auto *fsg1_48 = buffer.data(fsg1 + 48);
    const auto *fsg1_75 = buffer.data(fsg1 + 75);
    const auto *fsg1_80 = buffer.data(fsg1 + 80);
    const auto *fsg1_90 = buffer.data(fsg1 + 90);
    const auto *fsg1_93 = buffer.data(fsg1 + 93);
    const auto *fsg1_100 = buffer.data(fsg1 + 100);
    const auto *fsg1_102 = buffer.data(fsg1 + 102);
    const auto *fsg1_104 = buffer.data(fsg1 + 104);
    const auto *fsg1_110 = buffer.data(fsg1 + 110);
    const auto *fsg1_115 = buffer.data(fsg1 + 115);
    const auto *fsg1_117 = buffer.data(fsg1 + 117);
    const auto *fsg1_119 = buffer.data(fsg1 + 119);
    const auto *fsg1_123 = buffer.data(fsg1 + 123);
    const auto *fsg1_130 = buffer.data(fsg1 + 130);
    const auto *fsg1_132 = buffer.data(fsg1 + 132);

    const auto *gsd0_0 = buffer.data(gsd0 + 0);
    const auto *gsd0_3 = buffer.data(gsd0 + 3);
    const auto *gsd0_5 = buffer.data(gsd0 + 5);
    const auto *gsd0_9 = buffer.data(gsd0 + 9);
    const auto *gsd0_16 = buffer.data(gsd0 + 16);
    const auto *gsd0_17 = buffer.data(gsd0 + 17);
    const auto *gsd0_18 = buffer.data(gsd0 + 18);
    const auto *gsd0_21 = buffer.data(gsd0 + 21);
    const auto *gsd0_23 = buffer.data(gsd0 + 23);
    const auto *gsd0_29 = buffer.data(gsd0 + 29);
    const auto *gsd0_30 = buffer.data(gsd0 + 30);
    const auto *gsd0_33 = buffer.data(gsd0 + 33);
    const auto *gsd0_34 = buffer.data(gsd0 + 34);
    const auto *gsd0_35 = buffer.data(gsd0 + 35);
    const auto *gsd0_36 = buffer.data(gsd0 + 36);

    const auto *gsd1_0 = buffer.data(gsd1 + 0);
    const auto *gsd1_3 = buffer.data(gsd1 + 3);
    const auto *gsd1_5 = buffer.data(gsd1 + 5);
    const auto *gsd1_9 = buffer.data(gsd1 + 9);
    const auto *gsd1_16 = buffer.data(gsd1 + 16);
    const auto *gsd1_17 = buffer.data(gsd1 + 17);
    const auto *gsd1_18 = buffer.data(gsd1 + 18);
    const auto *gsd1_21 = buffer.data(gsd1 + 21);
    const auto *gsd1_23 = buffer.data(gsd1 + 23);
    const auto *gsd1_29 = buffer.data(gsd1 + 29);
    const auto *gsd1_30 = buffer.data(gsd1 + 30);
    const auto *gsd1_33 = buffer.data(gsd1 + 33);
    const auto *gsd1_34 = buffer.data(gsd1 + 34);
    const auto *gsd1_35 = buffer.data(gsd1 + 35);
    const auto *gsd1_36 = buffer.data(gsd1 + 36);

    const auto *gsf_0 = buffer.data(gsf + 0);
    const auto *gsf_1 = buffer.data(gsf + 1);
    const auto *gsf_2 = buffer.data(gsf + 2);
    const auto *gsf_3 = buffer.data(gsf + 3);
    const auto *gsf_5 = buffer.data(gsf + 5);
    const auto *gsf_6 = buffer.data(gsf + 6);
    const auto *gsf_8 = buffer.data(gsf + 8);
    const auto *gsf_9 = buffer.data(gsf + 9);
    const auto *gsf_10 = buffer.data(gsf + 10);
    const auto *gsf_11 = buffer.data(gsf + 11);
    const auto *gsf_13 = buffer.data(gsf + 13);
    const auto *gsf_16 = buffer.data(gsf + 16);
    const auto *gsf_17 = buffer.data(gsf + 17);
    const auto *gsf_18 = buffer.data(gsf + 18);
    const auto *gsf_19 = buffer.data(gsf + 19);
    const auto *gsf_20 = buffer.data(gsf + 20);
    const auto *gsf_22 = buffer.data(gsf + 22);
    const auto *gsf_25 = buffer.data(gsf + 25);
    const auto *gsf_26 = buffer.data(gsf + 26);
    const auto *gsf_27 = buffer.data(gsf + 27);
    const auto *gsf_28 = buffer.data(gsf + 28);
    const auto *gsf_29 = buffer.data(gsf + 29);
    const auto *gsf_30 = buffer.data(gsf + 30);
    const auto *gsf_31 = buffer.data(gsf + 31);
    const auto *gsf_32 = buffer.data(gsf + 32);
    const auto *gsf_33 = buffer.data(gsf + 33);
    const auto *gsf_36 = buffer.data(gsf + 36);
    const auto *gsf_37 = buffer.data(gsf + 37);
    const auto *gsf_38 = buffer.data(gsf + 38);
    const auto *gsf_39 = buffer.data(gsf + 39);
    const auto *gsf_40 = buffer.data(gsf + 40);
    const auto *gsf_42 = buffer.data(gsf + 42);
    const auto *gsf_46 = buffer.data(gsf + 46);
    const auto *gsf_47 = buffer.data(gsf + 47);
    const auto *gsf_48 = buffer.data(gsf + 48);
    const auto *gsf_49 = buffer.data(gsf + 49);
    const auto *gsf_50 = buffer.data(gsf + 50);
    const auto *gsf_51 = buffer.data(gsf + 51);
    const auto *gsf_52 = buffer.data(gsf + 52);
    const auto *gsf_55 = buffer.data(gsf + 55);
    const auto *gsf_56 = buffer.data(gsf + 56);
    const auto *gsf_57 = buffer.data(gsf + 57);
    const auto *gsf_58 = buffer.data(gsf + 58);
    const auto *gsf_59 = buffer.data(gsf + 59);
    const auto *gsf_60 = buffer.data(gsf + 60);
    const auto *gsf_61 = buffer.data(gsf + 61);
    const auto *gsf_62 = buffer.data(gsf + 62);
    const auto *gsf_63 = buffer.data(gsf + 63);
    const auto *gsf_66 = buffer.data(gsf + 66);
    const auto *gsf_68 = buffer.data(gsf + 68);
    const auto *gsf_69 = buffer.data(gsf + 69);
    const auto *gsf_70 = buffer.data(gsf + 70);
    const auto *gsf_72 = buffer.data(gsf + 72);
    const auto *gsf_76 = buffer.data(gsf + 76);
    const auto *gsf_77 = buffer.data(gsf + 77);
    const auto *gsf_78 = buffer.data(gsf + 78);
    const auto *gsf_79 = buffer.data(gsf + 79);
    const auto *gsf_80 = buffer.data(gsf + 80);
    const auto *gsf_82 = buffer.data(gsf + 82);
    const auto *gsf_86 = buffer.data(gsf + 86);
    const auto *gsf_87 = buffer.data(gsf + 87);
    const auto *gsf_88 = buffer.data(gsf + 88);
    const auto *gsf_89 = buffer.data(gsf + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, fsf_0, gsd0_0, \
                         gsd1_0, gsf_0, gsf_1, gsf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fsf_0[k]
                 + f_1 * gsd0_0[k]
                 - f_2 * gsd1_0[k]
                 + f_3 * pc_x[k] * gsf_0[k];

        t_1[k] = f_3 * pc_y[k] * gsf_0[k];

        t_2[k] = f_3 * pc_z[k] * gsf_0[k];

        t_3[k] = f_4 * gsd0_0[k]
                 - f_5 * gsd1_0[k]
                 + f_3 * pc_y[k] * gsf_1[k];

        t_4[k] = f_3 * pc_y[k] * gsf_2[k];

        t_5[k] = f_4 * gsd0_0[k]
                 - f_5 * gsd1_0[k]
                 + f_3 * pc_z[k] * gsf_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_x, pc_y, pc_z, fsf_6, fsf_9, gsd0_3, \
                         gsd1_3, gsf_3, gsf_5, gsf_6, gsf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * fsf_6[k]
                 + f_3 * pc_x[k] * gsf_6[k];

        t_7[k] = f_3 * pc_z[k] * gsf_3[k];

        t_8[k] = f_3 * pc_y[k] * gsf_5[k];

        t_9[k] = f_0 * fsf_9[k]
                 + f_3 * pc_x[k] * gsf_9[k];

        t_10[k] = f_1 * gsd0_3[k]
                  - f_2 * gsd1_3[k]
                  + f_3 * pc_y[k] * gsf_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pc_y, pc_z, fsg0_0, fsg1_0, \
                         gsd0_5, gsd1_5, gsf_6, gsf_8, gsf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * gsf_6[k];

        t_12[k] = f_4 * gsd0_5[k]
                  - f_5 * gsd1_5[k]
                  + f_3 * pc_y[k] * gsf_8[k];

        t_13[k] = f_3 * pc_y[k] * gsf_9[k];

        t_14[k] = f_1 * gsd0_5[k]
                  - f_2 * gsd1_5[k]
                  + f_3 * pc_z[k] * gsf_9[k];

        t_15[k] = pa_y[k] * fsg0_0[k]
                  - f_6 * pc_y[k] * fsg1_0[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_y, pc_y, pc_z, fsg0_3, fsg0_5, \
                         fsf_0, fsf_1, fsg1_3, fsg1_5, gsf_10, gsf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_7 * fsf_0[k]
                  + f_3 * pc_y[k] * gsf_10[k];

        t_17[k] = f_3 * pc_z[k] * gsf_10[k];

        t_18[k] = pa_y[k] * fsg0_3[k]
                  + f_8 * fsf_1[k]
                  - f_6 * pc_y[k] * fsg1_3[k];

        t_19[k] = f_3 * pc_z[k] * gsf_11[k];

        t_20[k] = pa_y[k] * fsg0_5[k]
                  - f_6 * pc_y[k] * fsg1_5[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_x, pc_z, fsf_16, fsf_18, fsf_19, gsf_13, \
                         gsf_16, gsf_18, gsf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_9 * fsf_16[k]
                  + f_3 * pc_x[k] * gsf_16[k];

        t_22[k] = f_3 * pc_z[k] * gsf_13[k];

        t_23[k] = f_9 * fsf_18[k]
                  + f_3 * pc_x[k] * gsf_18[k];

        t_24[k] = f_9 * fsf_19[k]
                  + f_3 * pc_x[k] * gsf_19[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pc_y, pc_z, fsf_6, fsf_9, gsd0_9, gsd1_9, \
                         gsf_16, gsf_17, gsf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_7 * fsf_6[k]
                  + f_1 * gsd0_9[k]
                  - f_2 * gsd1_9[k]
                  + f_3 * pc_y[k] * gsf_16[k];

        t_26[k] = f_3 * pc_z[k] * gsf_16[k];

        t_27[k] = f_4 * gsd0_9[k]
                  - f_5 * gsd1_9[k]
                  + f_3 * pc_z[k] * gsf_17[k];

        t_28[k] = f_7 * fsf_9[k]
                  + f_3 * pc_y[k] * gsf_19[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_y, pa_z, pc_y, pc_z, fsg0_0, fsg0_14, \
                         fsf_0, fsg1_0, fsg1_14, gsf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_y[k] * fsg0_14[k]
                  - f_6 * pc_y[k] * fsg1_14[k];

        t_30[k] = pa_z[k] * fsg0_0[k]
                  - f_6 * pc_z[k] * fsg1_0[k];

        t_31[k] = f_3 * pc_y[k] * gsf_20[k];

        t_32[k] = f_7 * fsf_0[k]
                  + f_3 * pc_z[k] * gsf_20[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_z, pc_x, pc_y, pc_z, fsg0_3, fsg0_5, \
                         fsf_2, fsf_26, fsg1_3, fsg1_5, gsf_22, \
                         gsf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_z[k] * fsg0_3[k]
                  - f_6 * pc_z[k] * fsg1_3[k];

        t_34[k] = f_3 * pc_y[k] * gsf_22[k];

        t_35[k] = pa_z[k] * fsg0_5[k]
                  + f_8 * fsf_2[k]
                  - f_6 * pc_z[k] * fsg1_5[k];

        t_36[k] = f_9 * fsf_26[k]
                  + f_3 * pc_x[k] * gsf_26[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_z, pc_x, pc_y, pc_z, fsg0_10, fsf_27, \
                         fsf_29, fsg1_10, gsf_25, gsf_27, gsf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_9 * fsf_27[k]
                  + f_3 * pc_x[k] * gsf_27[k];

        t_38[k] = f_3 * pc_y[k] * gsf_25[k];

        t_39[k] = f_9 * fsf_29[k]
                  + f_3 * pc_x[k] * gsf_29[k];

        t_40[k] = pa_z[k] * fsg0_10[k]
                  - f_6 * pc_z[k] * fsg1_10[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pc_y, pc_z, fsf_9, gsd0_16, gsd0_17, gsd1_16, \
                         gsd1_17, gsf_27, gsf_28, gsf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_10 * gsd0_16[k]
                  - f_11 * gsd1_16[k]
                  + f_3 * pc_y[k] * gsf_27[k];

        t_42[k] = f_4 * gsd0_17[k]
                  - f_5 * gsd1_17[k]
                  + f_3 * pc_y[k] * gsf_28[k];

        t_43[k] = f_3 * pc_y[k] * gsf_29[k];

        t_44[k] = f_7 * fsf_9[k]
                  + f_1 * gsd0_17[k]
                  - f_2 * gsd1_17[k]
                  + f_3 * pc_z[k] * gsf_29[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pc_x, pc_y, pc_z, fsf_10, fsf_30, fsf_33, \
                         gsd0_18, gsd0_21, gsd1_18, gsd1_21, gsf_30, \
                         gsf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_8 * fsf_30[k]
                  + f_1 * gsd0_18[k]
                  - f_2 * gsd1_18[k]
                  + f_3 * pc_x[k] * gsf_30[k];

        t_46[k] = f_8 * fsf_10[k]
                  + f_3 * pc_y[k] * gsf_30[k];

        t_47[k] = f_3 * pc_z[k] * gsf_30[k];

        t_48[k] = f_8 * fsf_33[k]
                  + f_4 * gsd0_21[k]
                  - f_5 * gsd1_21[k]
                  + f_3 * pc_x[k] * gsf_33[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, pc_x, pc_z, fsf_36, fsf_38, gsd0_18, \
                         gsd1_18, gsf_31, gsf_32, gsf_33, gsf_36, \
                         gsf_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_3 * pc_z[k] * gsf_31[k];

        t_50[k] = f_4 * gsd0_18[k]
                  - f_5 * gsd1_18[k]
                  + f_3 * pc_z[k] * gsf_32[k];

        t_51[k] = f_8 * fsf_36[k]
                  + f_3 * pc_x[k] * gsf_36[k];

        t_52[k] = f_3 * pc_z[k] * gsf_33[k];

        t_53[k] = f_8 * fsf_38[k]
                  + f_3 * pc_x[k] * gsf_38[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pc_x, pc_y, pc_z, fsf_16, fsf_19, \
                         fsf_39, gsd0_21, gsd1_21, gsf_36, gsf_37, \
                         gsf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_8 * fsf_39[k]
                  + f_3 * pc_x[k] * gsf_39[k];

        t_55[k] = f_8 * fsf_16[k]
                  + f_1 * gsd0_21[k]
                  - f_2 * gsd1_21[k]
                  + f_3 * pc_y[k] * gsf_36[k];

        t_56[k] = f_3 * pc_z[k] * gsf_36[k];

        t_57[k] = f_4 * gsd0_21[k]
                  - f_5 * gsd1_21[k]
                  + f_3 * pc_z[k] * gsf_37[k];

        t_58[k] = f_8 * fsf_19[k]
                  + f_3 * pc_y[k] * gsf_39[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_y, pc_y, pc_z, fsg0_30, fsf_10, fsf_20, \
                         fsg1_30, gsd0_23, gsd1_23, gsf_39, gsf_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_1 * gsd0_23[k]
                  - f_2 * gsd1_23[k]
                  + f_3 * pc_z[k] * gsf_39[k];

        t_60[k] = pa_y[k] * fsg0_30[k]
                  - f_6 * pc_y[k] * fsg1_30[k];

        t_61[k] = f_7 * fsf_20[k]
                  + f_3 * pc_y[k] * gsf_40[k];

        t_62[k] = f_7 * fsf_10[k]
                  + f_3 * pc_z[k] * gsf_40[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pa_y, pa_z, pc_y, pc_z, fsg0_18, fsg0_35, fsf_22, \
                         fsg1_18, fsg1_35, gsf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pa_z[k] * fsg0_18[k]
                  - f_6 * pc_z[k] * fsg1_18[k];

        t_64[k] = f_7 * fsf_22[k]
                  + f_3 * pc_y[k] * gsf_42[k];

        t_65[k] = pa_y[k] * fsg0_35[k]
                  - f_6 * pc_y[k] * fsg1_35[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pc_x, fsf_46, fsf_47, fsf_48, fsf_49, gsf_46, \
                         gsf_47, gsf_48, gsf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_8 * fsf_46[k]
                  + f_3 * pc_x[k] * gsf_46[k];

        t_67[k] = f_8 * fsf_47[k]
                  + f_3 * pc_x[k] * gsf_47[k];

        t_68[k] = f_8 * fsf_48[k]
                  + f_3 * pc_x[k] * gsf_48[k];

        t_69[k] = f_8 * fsf_49[k]
                  + f_3 * pc_x[k] * gsf_49[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pa_z, pc_y, pc_z, fsg0_25, fsf_16, fsf_28, fsg1_25, \
                         gsd0_29, gsd1_29, gsf_46, gsf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pa_z[k] * fsg0_25[k]
                  - f_6 * pc_z[k] * fsg1_25[k];

        t_71[k] = f_7 * fsf_16[k]
                  + f_3 * pc_z[k] * gsf_46[k];

        t_72[k] = f_7 * fsf_28[k]
                  + f_4 * gsd0_29[k]
                  - f_5 * gsd1_29[k]
                  + f_3 * pc_y[k] * gsf_48[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_y, pc_x, pc_y, fsg0_44, fsf_29, fsf_50, \
                         fsg1_44, gsd0_30, gsd1_30, gsf_49, gsf_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_7 * fsf_29[k]
                  + f_3 * pc_y[k] * gsf_49[k];

        t_74[k] = pa_y[k] * fsg0_44[k]
                  - f_6 * pc_y[k] * fsg1_44[k];

        t_75[k] = f_8 * fsf_50[k]
                  + f_1 * gsd0_30[k]
                  - f_2 * gsd1_30[k]
                  + f_3 * pc_x[k] * gsf_50[k];

        t_76[k] = f_3 * pc_y[k] * gsf_50[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pc_y, pc_z, fsf_20, gsd0_30, gsd1_30, gsf_50, \
                         gsf_51, gsf_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_8 * fsf_20[k]
                  + f_3 * pc_z[k] * gsf_50[k];

        t_78[k] = f_4 * gsd0_30[k]
                  - f_5 * gsd1_30[k]
                  + f_3 * pc_y[k] * gsf_51[k];

        t_79[k] = f_3 * pc_y[k] * gsf_52[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pc_x, pc_y, fsf_55, fsf_56, fsf_57, gsd0_35, \
                         gsd1_35, gsf_55, gsf_56, gsf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_8 * fsf_55[k]
                  + f_4 * gsd0_35[k]
                  - f_5 * gsd1_35[k]
                  + f_3 * pc_x[k] * gsf_55[k];

        t_81[k] = f_8 * fsf_56[k]
                  + f_3 * pc_x[k] * gsf_56[k];

        t_82[k] = f_8 * fsf_57[k]
                  + f_3 * pc_x[k] * gsf_57[k];

        t_83[k] = f_3 * pc_y[k] * gsf_55[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pc_x, pc_y, fsf_59, gsd0_33, gsd0_34, gsd1_33, \
                         gsd1_34, gsf_56, gsf_57, gsf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_8 * fsf_59[k]
                  + f_3 * pc_x[k] * gsf_59[k];

        t_85[k] = f_1 * gsd0_33[k]
                  - f_2 * gsd1_33[k]
                  + f_3 * pc_y[k] * gsf_56[k];

        t_86[k] = f_10 * gsd0_34[k]
                  - f_11 * gsd1_34[k]
                  + f_3 * pc_y[k] * gsf_57[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_x, pc_x, pc_y, pc_z, fsg0_90, fsf_29, \
                         fsf_60, fsg1_90, gsd0_35, gsd1_35, gsf_58, \
                         gsf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_4 * gsd0_35[k]
                  - f_5 * gsd1_35[k]
                  + f_3 * pc_y[k] * gsf_58[k];

        t_88[k] = f_3 * pc_y[k] * gsf_59[k];

        t_89[k] = f_8 * fsf_29[k]
                  + f_1 * gsd0_35[k]
                  - f_2 * gsd1_35[k]
                  + f_3 * pc_z[k] * gsf_59[k];

        t_90[k] = pa_x[k] * fsg0_90[k]
                  + f_0 * fsf_60[k]
                  - f_6 * pc_x[k] * fsg1_90[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pa_x, pc_x, pc_y, pc_z, fsg0_93, fsf_30, \
                         fsf_63, fsg1_93, gsf_60, gsf_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_9 * fsf_30[k]
                  + f_3 * pc_y[k] * gsf_60[k];

        t_92[k] = f_3 * pc_z[k] * gsf_60[k];

        t_93[k] = pa_x[k] * fsg0_93[k]
                  + f_8 * fsf_63[k]
                  - f_6 * pc_x[k] * fsg1_93[k];

        t_94[k] = f_3 * pc_z[k] * gsf_61[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pc_x, pc_z, fsf_66, fsf_68, gsd0_36, gsd1_36, \
                         gsf_62, gsf_63, gsf_66, gsf_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_4 * gsd0_36[k]
                  - f_5 * gsd1_36[k]
                  + f_3 * pc_z[k] * gsf_62[k];

        t_96[k] = f_7 * fsf_66[k]
                  + f_3 * pc_x[k] * gsf_66[k];

        t_97[k] = f_3 * pc_z[k] * gsf_63[k];

        t_98[k] = f_7 * fsf_68[k]
                  + f_3 * pc_x[k] * gsf_68[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pa_x, pc_x, pc_z, fsg0_100, fsg0_102, \
                         fsf_69, fsg1_100, fsg1_102, gsf_66, gsf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_7 * fsf_69[k]
                  + f_3 * pc_x[k] * gsf_69[k];

        t_100[k] = pa_x[k] * fsg0_100[k]
                   - f_6 * pc_x[k] * fsg1_100[k];

        t_101[k] = f_3 * pc_z[k] * gsf_66[k];

        t_102[k] = pa_x[k] * fsg0_102[k]
                   - f_6 * pc_x[k] * fsg1_102[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, pa_x, pa_z, pc_x, pc_y, pc_z, fsg0_45, fsg0_104, \
                         fsf_39, fsg1_45, fsg1_104, gsf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_9 * fsf_39[k]
                   + f_3 * pc_y[k] * gsf_69[k];

        t_104[k] = pa_x[k] * fsg0_104[k]
                   - f_6 * pc_x[k] * fsg1_104[k];

        t_105[k] = pa_z[k] * fsg0_45[k]
                   - f_6 * pc_z[k] * fsg1_45[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pa_z, pc_y, pc_z, fsg0_48, fsf_30, \
                         fsf_40, fsf_42, fsg1_48, gsf_70, gsf_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_8 * fsf_40[k]
                   + f_3 * pc_y[k] * gsf_70[k];

        t_107[k] = f_7 * fsf_30[k]
                   + f_3 * pc_z[k] * gsf_70[k];

        t_108[k] = pa_z[k] * fsg0_48[k]
                   - f_6 * pc_z[k] * fsg1_48[k];

        t_109[k] = f_8 * fsf_42[k]
                   + f_3 * pc_y[k] * gsf_72[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pa_x, pc_x, fsg0_110, fsf_75, fsf_76, \
                         fsf_77, fsf_78, fsg1_110, gsf_76, gsf_77, \
                         gsf_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pa_x[k] * fsg0_110[k]
                   + f_8 * fsf_75[k]
                   - f_6 * pc_x[k] * fsg1_110[k];

        t_111[k] = f_7 * fsf_76[k]
                   + f_3 * pc_x[k] * gsf_76[k];

        t_112[k] = f_7 * fsf_77[k]
                   + f_3 * pc_x[k] * gsf_77[k];

        t_113[k] = f_7 * fsf_78[k]
                   + f_3 * pc_x[k] * gsf_78[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pa_x, pc_x, pc_z, fsg0_115, fsg0_117, \
                         fsf_36, fsf_79, fsg1_115, fsg1_117, gsf_76, \
                         gsf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_7 * fsf_79[k]
                   + f_3 * pc_x[k] * gsf_79[k];

        t_115[k] = pa_x[k] * fsg0_115[k]
                   - f_6 * pc_x[k] * fsg1_115[k];

        t_116[k] = f_7 * fsf_36[k]
                   + f_3 * pc_z[k] * gsf_76[k];

        t_117[k] = pa_x[k] * fsg0_117[k]
                   - f_6 * pc_x[k] * fsg1_117[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pa_x, pa_y, pc_x, pc_y, fsg0_75, \
                         fsg0_119, fsf_49, fsf_50, fsg1_75, fsg1_119, gsf_79, \
                         gsf_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_8 * fsf_49[k]
                   + f_3 * pc_y[k] * gsf_79[k];

        t_119[k] = pa_x[k] * fsg0_119[k]
                   - f_6 * pc_x[k] * fsg1_119[k];

        t_120[k] = pa_y[k] * fsg0_75[k]
                   - f_6 * pc_y[k] * fsg1_75[k];

        t_121[k] = f_7 * fsf_50[k]
                   + f_3 * pc_y[k] * gsf_80[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, pa_x, pc_x, pc_y, pc_z, fsg0_123, fsf_40, \
                         fsf_52, fsf_83, fsg1_123, gsf_80, gsf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_8 * fsf_40[k]
                   + f_3 * pc_z[k] * gsf_80[k];

        t_123[k] = pa_x[k] * fsg0_123[k]
                   + f_8 * fsf_83[k]
                   - f_6 * pc_x[k] * fsg1_123[k];

        t_124[k] = f_7 * fsf_52[k]
                   + f_3 * pc_y[k] * gsf_82[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pa_y, pc_x, pc_y, fsg0_80, fsf_86, \
                         fsf_87, fsf_88, fsg1_80, gsf_86, gsf_87, \
                         gsf_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = pa_y[k] * fsg0_80[k]
                   - f_6 * pc_y[k] * fsg1_80[k];

        t_126[k] = f_7 * fsf_86[k]
                   + f_3 * pc_x[k] * gsf_86[k];

        t_127[k] = f_7 * fsf_87[k]
                   + f_3 * pc_x[k] * gsf_87[k];

        t_128[k] = f_7 * fsf_88[k]
                   + f_3 * pc_x[k] * gsf_88[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, pa_x, pc_x, pc_z, fsg0_130, fsg0_132, \
                         fsf_46, fsf_89, fsg1_130, fsg1_132, gsf_86, \
                         gsf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_7 * fsf_89[k]
                   + f_3 * pc_x[k] * gsf_89[k];

        t_130[k] = pa_x[k] * fsg0_130[k]
                   - f_6 * pc_x[k] * fsg1_130[k];

        t_131[k] = f_8 * fsf_46[k]
                   + f_3 * pc_z[k] * gsf_86[k];

        t_132[k] = pa_x[k] * fsg0_132[k]
                   - f_6 * pc_x[k] * fsg1_132[k];
    }
}

static auto
compute_prim_gsg_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t fsg0,
                                                          const size_t fsf, const size_t fsg1,
                                                          const size_t gsd0, const size_t gsd1,
                                                          const size_t gsf, const size_t ncols,
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
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *fsg0_90 = buffer.data(fsg0 + 90);
    const auto *fsg0_91 = buffer.data(fsg0 + 91);
    const auto *fsg0_93 = buffer.data(fsg0 + 93);
    const auto *fsg0_100 = buffer.data(fsg0 + 100);
    const auto *fsg0_102 = buffer.data(fsg0 + 102);
    const auto *fsg0_134 = buffer.data(fsg0 + 134);
    const auto *fsg0_135 = buffer.data(fsg0 + 135);
    const auto *fsg0_137 = buffer.data(fsg0 + 137);
    const auto *fsg0_140 = buffer.data(fsg0 + 140);
    const auto *fsg0_145 = buffer.data(fsg0 + 145);
    const auto *fsg0_146 = buffer.data(fsg0 + 146);
    const auto *fsg0_147 = buffer.data(fsg0 + 147);
    const auto *fsg0_149 = buffer.data(fsg0 + 149);

    const auto *fsf_50 = buffer.data(fsf + 50);
    const auto *fsf_59 = buffer.data(fsf + 59);
    const auto *fsf_66 = buffer.data(fsf + 66);
    const auto *fsf_67 = buffer.data(fsf + 67);
    const auto *fsf_69 = buffer.data(fsf + 69);
    const auto *fsf_76 = buffer.data(fsf + 76);
    const auto *fsf_79 = buffer.data(fsf + 79);
    const auto *fsf_86 = buffer.data(fsf + 86);
    const auto *fsf_88 = buffer.data(fsf + 88);
    const auto *fsf_89 = buffer.data(fsf + 89);
    const auto *fsf_90 = buffer.data(fsf + 90);
    const auto *fsf_95 = buffer.data(fsf + 95);
    const auto *fsf_96 = buffer.data(fsf + 96);
    const auto *fsf_97 = buffer.data(fsf + 97);
    const auto *fsf_98 = buffer.data(fsf + 98);
    const auto *fsf_99 = buffer.data(fsf + 99);

    const auto *fsg1_90 = buffer.data(fsg1 + 90);
    const auto *fsg1_91 = buffer.data(fsg1 + 91);
    const auto *fsg1_93 = buffer.data(fsg1 + 93);
    const auto *fsg1_100 = buffer.data(fsg1 + 100);
    const auto *fsg1_102 = buffer.data(fsg1 + 102);
    const auto *fsg1_134 = buffer.data(fsg1 + 134);
    const auto *fsg1_135 = buffer.data(fsg1 + 135);
    const auto *fsg1_137 = buffer.data(fsg1 + 137);
    const auto *fsg1_140 = buffer.data(fsg1 + 140);
    const auto *fsg1_145 = buffer.data(fsg1 + 145);
    const auto *fsg1_146 = buffer.data(fsg1 + 146);
    const auto *fsg1_147 = buffer.data(fsg1 + 147);
    const auto *fsg1_149 = buffer.data(fsg1 + 149);

    const auto *gsd0_54 = buffer.data(gsd0 + 54);
    const auto *gsd0_60 = buffer.data(gsd0 + 60);
    const auto *gsd0_61 = buffer.data(gsd0 + 61);
    const auto *gsd0_63 = buffer.data(gsd0 + 63);
    const auto *gsd0_65 = buffer.data(gsd0 + 65);
    const auto *gsd0_68 = buffer.data(gsd0 + 68);
    const auto *gsd0_70 = buffer.data(gsd0 + 70);
    const auto *gsd0_71 = buffer.data(gsd0 + 71);
    const auto *gsd0_72 = buffer.data(gsd0 + 72);
    const auto *gsd0_73 = buffer.data(gsd0 + 73);
    const auto *gsd0_74 = buffer.data(gsd0 + 74);
    const auto *gsd0_75 = buffer.data(gsd0 + 75);
    const auto *gsd0_76 = buffer.data(gsd0 + 76);
    const auto *gsd0_77 = buffer.data(gsd0 + 77);
    const auto *gsd0_79 = buffer.data(gsd0 + 79);
    const auto *gsd0_81 = buffer.data(gsd0 + 81);
    const auto *gsd0_82 = buffer.data(gsd0 + 82);
    const auto *gsd0_84 = buffer.data(gsd0 + 84);
    const auto *gsd0_86 = buffer.data(gsd0 + 86);
    const auto *gsd0_87 = buffer.data(gsd0 + 87);
    const auto *gsd0_88 = buffer.data(gsd0 + 88);
    const auto *gsd0_89 = buffer.data(gsd0 + 89);

    const auto *gsd1_54 = buffer.data(gsd1 + 54);
    const auto *gsd1_60 = buffer.data(gsd1 + 60);
    const auto *gsd1_61 = buffer.data(gsd1 + 61);
    const auto *gsd1_63 = buffer.data(gsd1 + 63);
    const auto *gsd1_65 = buffer.data(gsd1 + 65);
    const auto *gsd1_68 = buffer.data(gsd1 + 68);
    const auto *gsd1_70 = buffer.data(gsd1 + 70);
    const auto *gsd1_71 = buffer.data(gsd1 + 71);
    const auto *gsd1_72 = buffer.data(gsd1 + 72);
    const auto *gsd1_73 = buffer.data(gsd1 + 73);
    const auto *gsd1_74 = buffer.data(gsd1 + 74);
    const auto *gsd1_75 = buffer.data(gsd1 + 75);
    const auto *gsd1_76 = buffer.data(gsd1 + 76);
    const auto *gsd1_77 = buffer.data(gsd1 + 77);
    const auto *gsd1_79 = buffer.data(gsd1 + 79);
    const auto *gsd1_81 = buffer.data(gsd1 + 81);
    const auto *gsd1_82 = buffer.data(gsd1 + 82);
    const auto *gsd1_84 = buffer.data(gsd1 + 84);
    const auto *gsd1_86 = buffer.data(gsd1 + 86);
    const auto *gsd1_87 = buffer.data(gsd1 + 87);
    const auto *gsd1_88 = buffer.data(gsd1 + 88);
    const auto *gsd1_89 = buffer.data(gsd1 + 89);

    const auto *gsf_89 = buffer.data(gsf + 89);
    const auto *gsf_90 = buffer.data(gsf + 90);
    const auto *gsf_91 = buffer.data(gsf + 91);
    const auto *gsf_92 = buffer.data(gsf + 92);
    const auto *gsf_95 = buffer.data(gsf + 95);
    const auto *gsf_96 = buffer.data(gsf + 96);
    const auto *gsf_97 = buffer.data(gsf + 97);
    const auto *gsf_99 = buffer.data(gsf + 99);
    const auto *gsf_100 = buffer.data(gsf + 100);
    const auto *gsf_101 = buffer.data(gsf + 101);
    const auto *gsf_103 = buffer.data(gsf + 103);
    const auto *gsf_105 = buffer.data(gsf + 105);
    const auto *gsf_106 = buffer.data(gsf + 106);
    const auto *gsf_107 = buffer.data(gsf + 107);
    const auto *gsf_108 = buffer.data(gsf + 108);
    const auto *gsf_109 = buffer.data(gsf + 109);
    const auto *gsf_112 = buffer.data(gsf + 112);
    const auto *gsf_114 = buffer.data(gsf + 114);
    const auto *gsf_115 = buffer.data(gsf + 115);
    const auto *gsf_116 = buffer.data(gsf + 116);
    const auto *gsf_117 = buffer.data(gsf + 117);
    const auto *gsf_118 = buffer.data(gsf + 118);
    const auto *gsf_119 = buffer.data(gsf + 119);
    const auto *gsf_120 = buffer.data(gsf + 120);
    const auto *gsf_121 = buffer.data(gsf + 121);
    const auto *gsf_122 = buffer.data(gsf + 122);
    const auto *gsf_123 = buffer.data(gsf + 123);
    const auto *gsf_124 = buffer.data(gsf + 124);
    const auto *gsf_125 = buffer.data(gsf + 125);
    const auto *gsf_126 = buffer.data(gsf + 126);
    const auto *gsf_127 = buffer.data(gsf + 127);
    const auto *gsf_128 = buffer.data(gsf + 128);
    const auto *gsf_129 = buffer.data(gsf + 129);
    const auto *gsf_131 = buffer.data(gsf + 131);
    const auto *gsf_133 = buffer.data(gsf + 133);
    const auto *gsf_134 = buffer.data(gsf + 134);
    const auto *gsf_136 = buffer.data(gsf + 136);
    const auto *gsf_137 = buffer.data(gsf + 137);
    const auto *gsf_138 = buffer.data(gsf + 138);
    const auto *gsf_139 = buffer.data(gsf + 139);
    const auto *gsf_140 = buffer.data(gsf + 140);
    const auto *gsf_142 = buffer.data(gsf + 142);
    const auto *gsf_143 = buffer.data(gsf + 143);
    const auto *gsf_145 = buffer.data(gsf + 145);
    const auto *gsf_146 = buffer.data(gsf + 146);
    const auto *gsf_147 = buffer.data(gsf + 147);
    const auto *gsf_148 = buffer.data(gsf + 148);
    const auto *gsf_149 = buffer.data(gsf + 149);

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pa_x, pc_x, pc_y, fsg0_134, fsg0_135, \
                         fsf_59, fsf_90, fsg1_134, fsg1_135, gsf_89, \
                         gsf_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_7 * fsf_59[k]
                   + f_3 * pc_y[k] * gsf_89[k];

        t_134[k] = pa_x[k] * fsg0_134[k]
                   - f_6 * pc_x[k] * fsg1_134[k];

        t_135[k] = pa_x[k] * fsg0_135[k]
                   + f_0 * fsf_90[k]
                   - f_6 * pc_x[k] * fsg1_135[k];

        t_136[k] = f_3 * pc_y[k] * gsf_90[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pc_y, pc_z, fsf_50, gsd0_54, gsd1_54, gsf_90, \
                         gsf_91, gsf_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_9 * fsf_50[k]
                   + f_3 * pc_z[k] * gsf_90[k];

        t_138[k] = f_4 * gsd0_54[k]
                   - f_5 * gsd1_54[k]
                   + f_3 * pc_y[k] * gsf_91[k];

        t_139[k] = f_3 * pc_y[k] * gsf_92[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pa_x, pc_x, pc_y, fsg0_140, fsf_95, \
                         fsf_96, fsf_97, fsg1_140, gsf_95, gsf_96, \
                         gsf_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = pa_x[k] * fsg0_140[k]
                   + f_8 * fsf_95[k]
                   - f_6 * pc_x[k] * fsg1_140[k];

        t_141[k] = f_7 * fsf_96[k]
                   + f_3 * pc_x[k] * gsf_96[k];

        t_142[k] = f_7 * fsf_97[k]
                   + f_3 * pc_x[k] * gsf_97[k];

        t_143[k] = f_3 * pc_y[k] * gsf_95[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, pa_x, pc_x, pc_y, fsg0_145, \
                         fsg0_146, fsg0_147, fsf_99, fsg1_145, fsg1_146, fsg1_147, \
                         gsf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_7 * fsf_99[k]
                   + f_3 * pc_x[k] * gsf_99[k];

        t_145[k] = pa_x[k] * fsg0_145[k]
                   - f_6 * pc_x[k] * fsg1_145[k];

        t_146[k] = pa_x[k] * fsg0_146[k]
                   - f_6 * pc_x[k] * fsg1_146[k];

        t_147[k] = pa_x[k] * fsg0_147[k]
                   - f_6 * pc_x[k] * fsg1_147[k];

        t_148[k] = f_3 * pc_y[k] * gsf_99[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, pa_x, pc_x, pc_z, fsg0_149, fsg1_149, \
                         gsd0_60, gsd0_61, gsd1_60, gsd1_61, gsf_100, \
                         gsf_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = pa_x[k] * fsg0_149[k]
                   - f_6 * pc_x[k] * fsg1_149[k];

        t_150[k] = f_1 * gsd0_60[k]
                   - f_2 * gsd1_60[k]
                   + f_3 * pc_x[k] * gsf_100[k];

        t_151[k] = f_10 * gsd0_61[k]
                   - f_11 * gsd1_61[k]
                   + f_3 * pc_x[k] * gsf_101[k];

        t_152[k] = f_3 * pc_z[k] * gsf_100[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, t_157, pc_x, pc_z, gsd0_63, gsd0_65, \
                         gsd1_63, gsd1_65, gsf_101, gsf_103, gsf_105, gsf_106, \
                         gsf_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_4 * gsd0_63[k]
                   - f_5 * gsd1_63[k]
                   + f_3 * pc_x[k] * gsf_103[k];

        t_154[k] = f_3 * pc_z[k] * gsf_101[k];

        t_155[k] = f_4 * gsd0_65[k]
                   - f_5 * gsd1_65[k]
                   + f_3 * pc_x[k] * gsf_105[k];

        t_156[k] = f_3 * pc_x[k] * gsf_106[k];

        t_157[k] = f_3 * pc_x[k] * gsf_107[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, t_162, pc_x, pc_y, pc_z, fsf_66, gsd0_63, \
                         gsd1_63, gsf_106, gsf_107, gsf_108, gsf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_3 * pc_x[k] * gsf_108[k];

        t_159[k] = f_3 * pc_x[k] * gsf_109[k];

        t_160[k] = f_0 * fsf_66[k]
                   + f_1 * gsd0_63[k]
                   - f_2 * gsd1_63[k]
                   + f_3 * pc_y[k] * gsf_106[k];

        t_161[k] = f_3 * pc_z[k] * gsf_106[k];

        t_162[k] = f_4 * gsd0_63[k]
                   - f_5 * gsd1_63[k]
                   + f_3 * pc_z[k] * gsf_107[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, pa_z, pc_y, pc_z, fsg0_90, fsg0_91, \
                         fsf_69, fsg1_90, fsg1_91, gsd0_65, gsd1_65, \
                         gsf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_0 * fsf_69[k]
                   + f_3 * pc_y[k] * gsf_109[k];

        t_164[k] = f_1 * gsd0_65[k]
                   - f_2 * gsd1_65[k]
                   + f_3 * pc_z[k] * gsf_109[k];

        t_165[k] = pa_z[k] * fsg0_90[k]
                   - f_6 * pc_z[k] * fsg1_90[k];

        t_166[k] = pa_z[k] * fsg0_91[k]
                   - f_6 * pc_z[k] * fsg1_91[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, pa_z, pc_x, pc_z, fsg0_93, fsg1_93, gsd0_68, \
                         gsd0_70, gsd1_68, gsd1_70, gsf_112, gsf_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_10 * gsd0_68[k]
                   - f_11 * gsd1_68[k]
                   + f_3 * pc_x[k] * gsf_112[k];

        t_168[k] = pa_z[k] * fsg0_93[k]
                   - f_6 * pc_z[k] * fsg1_93[k];

        t_169[k] = f_4 * gsd0_70[k]
                   - f_5 * gsd1_70[k]
                   + f_3 * pc_x[k] * gsf_114[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, pc_x, gsd0_71, gsd1_71, gsf_115, \
                         gsf_116, gsf_117, gsf_118, gsf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_4 * gsd0_71[k]
                   - f_5 * gsd1_71[k]
                   + f_3 * pc_x[k] * gsf_115[k];

        t_171[k] = f_3 * pc_x[k] * gsf_116[k];

        t_172[k] = f_3 * pc_x[k] * gsf_117[k];

        t_173[k] = f_3 * pc_x[k] * gsf_118[k];

        t_174[k] = f_3 * pc_x[k] * gsf_119[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, pa_z, pc_y, pc_z, fsg0_100, fsg0_102, \
                         fsf_66, fsf_67, fsf_79, fsg1_100, fsg1_102, gsf_116, \
                         gsf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = pa_z[k] * fsg0_100[k]
                   - f_6 * pc_z[k] * fsg1_100[k];

        t_176[k] = f_7 * fsf_66[k]
                   + f_3 * pc_z[k] * gsf_116[k];

        t_177[k] = pa_z[k] * fsg0_102[k]
                   + f_8 * fsf_67[k]
                   - f_6 * pc_z[k] * fsg1_102[k];

        t_178[k] = f_9 * fsf_79[k]
                   + f_3 * pc_y[k] * gsf_119[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, pc_x, pc_z, fsf_69, gsd0_71, gsd0_72, gsd0_73, \
                         gsd1_71, gsd1_72, gsd1_73, gsf_119, gsf_120, \
                         gsf_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_7 * fsf_69[k]
                   + f_1 * gsd0_71[k]
                   - f_2 * gsd1_71[k]
                   + f_3 * pc_z[k] * gsf_119[k];

        t_180[k] = f_1 * gsd0_72[k]
                   - f_2 * gsd1_72[k]
                   + f_3 * pc_x[k] * gsf_120[k];

        t_181[k] = f_10 * gsd0_73[k]
                   - f_11 * gsd1_73[k]
                   + f_3 * pc_x[k] * gsf_121[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, pc_x, gsd0_74, gsd0_75, gsd0_76, gsd1_74, \
                         gsd1_75, gsd1_76, gsf_122, gsf_123, gsf_124 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = f_10 * gsd0_74[k]
                   - f_11 * gsd1_74[k]
                   + f_3 * pc_x[k] * gsf_122[k];

        t_183[k] = f_4 * gsd0_75[k]
                   - f_5 * gsd1_75[k]
                   + f_3 * pc_x[k] * gsf_123[k];

        t_184[k] = f_4 * gsd0_76[k]
                   - f_5 * gsd1_76[k]
                   + f_3 * pc_x[k] * gsf_124[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, t_189, pc_x, gsd0_77, gsd1_77, gsf_125, \
                         gsf_126, gsf_127, gsf_128, gsf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_4 * gsd0_77[k]
                   - f_5 * gsd1_77[k]
                   + f_3 * pc_x[k] * gsf_125[k];

        t_186[k] = f_3 * pc_x[k] * gsf_126[k];

        t_187[k] = f_3 * pc_x[k] * gsf_127[k];

        t_188[k] = f_3 * pc_x[k] * gsf_128[k];

        t_189[k] = f_3 * pc_x[k] * gsf_129[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, pc_y, pc_z, fsf_76, fsf_86, fsf_88, gsd0_75, \
                         gsd0_77, gsd1_75, gsd1_77, gsf_126, gsf_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = f_8 * fsf_86[k]
                   + f_1 * gsd0_75[k]
                   - f_2 * gsd1_75[k]
                   + f_3 * pc_y[k] * gsf_126[k];

        t_191[k] = f_8 * fsf_76[k]
                   + f_3 * pc_z[k] * gsf_126[k];

        t_192[k] = f_8 * fsf_88[k]
                   + f_4 * gsd0_77[k]
                   - f_5 * gsd1_77[k]
                   + f_3 * pc_y[k] * gsf_128[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, pa_y, pc_y, pc_z, fsg0_135, fsf_79, fsf_89, \
                         fsg1_135, gsd0_77, gsd1_77, gsf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_8 * fsf_89[k]
                   + f_3 * pc_y[k] * gsf_129[k];

        t_194[k] = f_8 * fsf_79[k]
                   + f_1 * gsd0_77[k]
                   - f_2 * gsd1_77[k]
                   + f_3 * pc_z[k] * gsf_129[k];

        t_195[k] = pa_y[k] * fsg0_135[k]
                   - f_6 * pc_y[k] * fsg1_135[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, pa_y, pc_x, pc_y, fsg0_137, fsg1_137, gsd0_79, \
                         gsd0_81, gsd1_79, gsd1_81, gsf_131, gsf_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_10 * gsd0_79[k]
                   - f_11 * gsd1_79[k]
                   + f_3 * pc_x[k] * gsf_131[k];

        t_197[k] = pa_y[k] * fsg0_137[k]
                   - f_6 * pc_y[k] * fsg1_137[k];

        t_198[k] = f_4 * gsd0_81[k]
                   - f_5 * gsd1_81[k]
                   + f_3 * pc_x[k] * gsf_133[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, t_203, pa_y, pc_x, pc_y, fsg0_140, \
                         fsg1_140, gsd0_82, gsd1_82, gsf_134, gsf_136, gsf_137, \
                         gsf_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_4 * gsd0_82[k]
                   - f_5 * gsd1_82[k]
                   + f_3 * pc_x[k] * gsf_134[k];

        t_200[k] = pa_y[k] * fsg0_140[k]
                   - f_6 * pc_y[k] * fsg1_140[k];

        t_201[k] = f_3 * pc_x[k] * gsf_136[k];

        t_202[k] = f_3 * pc_x[k] * gsf_137[k];

        t_203[k] = f_3 * pc_x[k] * gsf_138[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, pa_y, pc_x, pc_y, pc_z, fsg0_145, fsf_86, \
                         fsf_96, fsg1_145, gsf_136, gsf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_3 * pc_x[k] * gsf_139[k];

        t_205[k] = pa_y[k] * fsg0_145[k]
                   + f_0 * fsf_96[k]
                   - f_6 * pc_y[k] * fsg1_145[k];

        t_206[k] = f_9 * fsf_86[k]
                   + f_3 * pc_z[k] * gsf_136[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pa_y, pc_y, fsg0_147, fsg0_149, fsf_98, fsf_99, \
                         fsg1_147, fsg1_149, gsf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = pa_y[k] * fsg0_147[k]
                   + f_8 * fsf_98[k]
                   - f_6 * pc_y[k] * fsg1_147[k];

        t_208[k] = f_7 * fsf_99[k]
                   + f_3 * pc_y[k] * gsf_139[k];

        t_209[k] = pa_y[k] * fsg0_149[k]
                   - f_6 * pc_y[k] * fsg1_149[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, pc_x, pc_y, gsd0_84, gsd0_86, \
                         gsd0_87, gsd1_84, gsd1_86, gsd1_87, gsf_140, gsf_142, \
                         gsf_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_1 * gsd0_84[k]
                   - f_2 * gsd1_84[k]
                   + f_3 * pc_x[k] * gsf_140[k];

        t_211[k] = f_3 * pc_y[k] * gsf_140[k];

        t_212[k] = f_10 * gsd0_86[k]
                   - f_11 * gsd1_86[k]
                   + f_3 * pc_x[k] * gsf_142[k];

        t_213[k] = f_4 * gsd0_87[k]
                   - f_5 * gsd1_87[k]
                   + f_3 * pc_x[k] * gsf_143[k];

        t_214[k] = f_3 * pc_y[k] * gsf_142[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, t_219, pc_x, gsd0_89, gsd1_89, gsf_145, \
                         gsf_146, gsf_147, gsf_148, gsf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_4 * gsd0_89[k]
                   - f_5 * gsd1_89[k]
                   + f_3 * pc_x[k] * gsf_145[k];

        t_216[k] = f_3 * pc_x[k] * gsf_146[k];

        t_217[k] = f_3 * pc_x[k] * gsf_147[k];

        t_218[k] = f_3 * pc_x[k] * gsf_148[k];

        t_219[k] = f_3 * pc_x[k] * gsf_149[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, pc_y, gsd0_87, gsd0_88, gsd0_89, gsd1_87, \
                         gsd1_88, gsd1_89, gsf_146, gsf_147, gsf_148, \
                         gsf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_1 * gsd0_87[k]
                   - f_2 * gsd1_87[k]
                   + f_3 * pc_y[k] * gsf_146[k];

        t_221[k] = f_10 * gsd0_88[k]
                   - f_11 * gsd1_88[k]
                   + f_3 * pc_y[k] * gsf_147[k];

        t_222[k] = f_4 * gsd0_89[k]
                   - f_5 * gsd1_89[k]
                   + f_3 * pc_y[k] * gsf_148[k];

        t_223[k] = f_3 * pc_y[k] * gsf_149[k];
    }

#pragma omp simd aligned(t_224, pc_z, fsf_99, gsd0_89, gsd1_89, \
                         gsf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = f_0 * fsf_99[k]
                   + f_1 * gsd0_89[k]
                   - f_2 * gsd1_89[k]
                   + f_3 * pc_z[k] * gsf_149[k];
    }
}

auto
compute_prim_gsg_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t fsg0, const size_t fsf,
                                                   const size_t fsg1, const size_t gsd0,
                                                   const size_t gsd1, const size_t gsf,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_gsg_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, fsg0, fsf,
                                                              fsg1, gsd0, gsd1, gsf, ncols,
                                                              gamma, p, q);

    compute_prim_gsg_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, fsg0, fsf,
                                                              fsg1, gsd0, gsd1, gsf, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
