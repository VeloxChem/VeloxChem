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


#include "SimdThreeCenterElectronRepulsionVrrRecGPG.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_gpg_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t dpg0, const size_t dpg1,
                                                          const size_t fpg0, const size_t fpf,
                                                          const size_t fpg1, const size_t gsg0,
                                                          const size_t gsf, const size_t gsg1,
                                                          const size_t gpd0, const size_t gpd1,
                                                          const size_t gpf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 1.5 / gamma;
    const auto f_3 = 1.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = gamma / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 1.5 / q;
    const auto f_10 = 1.0 / p;
    const auto f_11 = gamma / (p * q);
    const auto f_12 = 1.0 / gamma;
    const auto f_13 = p / (gamma * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dpg0_70 = buffer.data(dpg0 + 70);

    const auto *dpg1_70 = buffer.data(dpg1 + 70);

    const auto *fpg0_0 = buffer.data(fpg0 + 0);
    const auto *fpg0_3 = buffer.data(fpg0 + 3);
    const auto *fpg0_5 = buffer.data(fpg0 + 5);
    const auto *fpg0_6 = buffer.data(fpg0 + 6);
    const auto *fpg0_9 = buffer.data(fpg0 + 9);
    const auto *fpg0_10 = buffer.data(fpg0 + 10);
    const auto *fpg0_14 = buffer.data(fpg0 + 14);
    const auto *fpg0_15 = buffer.data(fpg0 + 15);
    const auto *fpg0_18 = buffer.data(fpg0 + 18);
    const auto *fpg0_25 = buffer.data(fpg0 + 25);
    const auto *fpg0_30 = buffer.data(fpg0 + 30);
    const auto *fpg0_35 = buffer.data(fpg0 + 35);
    const auto *fpg0_44 = buffer.data(fpg0 + 44);
    const auto *fpg0_70 = buffer.data(fpg0 + 70);

    const auto *fpf_0 = buffer.data(fpf + 0);
    const auto *fpf_1 = buffer.data(fpf + 1);
    const auto *fpf_2 = buffer.data(fpf + 2);
    const auto *fpf_6 = buffer.data(fpf + 6);
    const auto *fpf_9 = buffer.data(fpf + 9);
    const auto *fpf_10 = buffer.data(fpf + 10);
    const auto *fpf_16 = buffer.data(fpf + 16);
    const auto *fpf_19 = buffer.data(fpf + 19);
    const auto *fpf_20 = buffer.data(fpf + 20);
    const auto *fpf_26 = buffer.data(fpf + 26);
    const auto *fpf_29 = buffer.data(fpf + 29);
    const auto *fpf_36 = buffer.data(fpf + 36);
    const auto *fpf_38 = buffer.data(fpf + 38);
    const auto *fpf_40 = buffer.data(fpf + 40);
    const auto *fpf_43 = buffer.data(fpf + 43);
    const auto *fpf_46 = buffer.data(fpf + 46);
    const auto *fpf_48 = buffer.data(fpf + 48);
    const auto *fpf_49 = buffer.data(fpf + 49);
    const auto *fpf_56 = buffer.data(fpf + 56);
    const auto *fpf_58 = buffer.data(fpf + 58);
    const auto *fpf_59 = buffer.data(fpf + 59);
    const auto *fpf_67 = buffer.data(fpf + 67);
    const auto *fpf_69 = buffer.data(fpf + 69);
    const auto *fpf_76 = buffer.data(fpf + 76);
    const auto *fpf_77 = buffer.data(fpf + 77);
    const auto *fpf_79 = buffer.data(fpf + 79);
    const auto *fpf_80 = buffer.data(fpf + 80);
    const auto *fpf_85 = buffer.data(fpf + 85);
    const auto *fpf_86 = buffer.data(fpf + 86);
    const auto *fpf_87 = buffer.data(fpf + 87);

    const auto *fpg1_0 = buffer.data(fpg1 + 0);
    const auto *fpg1_3 = buffer.data(fpg1 + 3);
    const auto *fpg1_5 = buffer.data(fpg1 + 5);
    const auto *fpg1_6 = buffer.data(fpg1 + 6);
    const auto *fpg1_9 = buffer.data(fpg1 + 9);
    const auto *fpg1_10 = buffer.data(fpg1 + 10);
    const auto *fpg1_14 = buffer.data(fpg1 + 14);
    const auto *fpg1_15 = buffer.data(fpg1 + 15);
    const auto *fpg1_18 = buffer.data(fpg1 + 18);
    const auto *fpg1_25 = buffer.data(fpg1 + 25);
    const auto *fpg1_30 = buffer.data(fpg1 + 30);
    const auto *fpg1_35 = buffer.data(fpg1 + 35);
    const auto *fpg1_44 = buffer.data(fpg1 + 44);
    const auto *fpg1_70 = buffer.data(fpg1 + 70);

    const auto *gsg0_0 = buffer.data(gsg0 + 0);
    const auto *gsg0_3 = buffer.data(gsg0 + 3);
    const auto *gsg0_5 = buffer.data(gsg0 + 5);
    const auto *gsg0_10 = buffer.data(gsg0 + 10);
    const auto *gsg0_12 = buffer.data(gsg0 + 12);
    const auto *gsg0_14 = buffer.data(gsg0 + 14);
    const auto *gsg0_18 = buffer.data(gsg0 + 18);
    const auto *gsg0_25 = buffer.data(gsg0 + 25);
    const auto *gsg0_27 = buffer.data(gsg0 + 27);
    const auto *gsg0_35 = buffer.data(gsg0 + 35);
    const auto *gsg0_41 = buffer.data(gsg0 + 41);
    const auto *gsg0_42 = buffer.data(gsg0 + 42);
    const auto *gsg0_44 = buffer.data(gsg0 + 44);

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
    const auto *gsf_27 = buffer.data(gsf + 27);
    const auto *gsf_28 = buffer.data(gsf + 28);
    const auto *gsf_29 = buffer.data(gsf + 29);

    const auto *gsg1_0 = buffer.data(gsg1 + 0);
    const auto *gsg1_3 = buffer.data(gsg1 + 3);
    const auto *gsg1_5 = buffer.data(gsg1 + 5);
    const auto *gsg1_10 = buffer.data(gsg1 + 10);
    const auto *gsg1_12 = buffer.data(gsg1 + 12);
    const auto *gsg1_14 = buffer.data(gsg1 + 14);
    const auto *gsg1_18 = buffer.data(gsg1 + 18);
    const auto *gsg1_25 = buffer.data(gsg1 + 25);
    const auto *gsg1_27 = buffer.data(gsg1 + 27);
    const auto *gsg1_35 = buffer.data(gsg1 + 35);
    const auto *gsg1_41 = buffer.data(gsg1 + 41);
    const auto *gsg1_42 = buffer.data(gsg1 + 42);
    const auto *gsg1_44 = buffer.data(gsg1 + 44);

    const auto *gpd0_0 = buffer.data(gpd0 + 0);
    const auto *gpd0_3 = buffer.data(gpd0 + 3);
    const auto *gpd0_5 = buffer.data(gpd0 + 5);
    const auto *gpd0_17 = buffer.data(gpd0 + 17);
    const auto *gpd0_21 = buffer.data(gpd0 + 21);
    const auto *gpd0_24 = buffer.data(gpd0 + 24);
    const auto *gpd0_27 = buffer.data(gpd0 + 27);
    const auto *gpd0_29 = buffer.data(gpd0 + 29);
    const auto *gpd0_40 = buffer.data(gpd0 + 40);
    const auto *gpd0_41 = buffer.data(gpd0 + 41);
    const auto *gpd0_48 = buffer.data(gpd0 + 48);
    const auto *gpd0_53 = buffer.data(gpd0 + 53);

    const auto *gpd1_0 = buffer.data(gpd1 + 0);
    const auto *gpd1_3 = buffer.data(gpd1 + 3);
    const auto *gpd1_5 = buffer.data(gpd1 + 5);
    const auto *gpd1_17 = buffer.data(gpd1 + 17);
    const auto *gpd1_21 = buffer.data(gpd1 + 21);
    const auto *gpd1_24 = buffer.data(gpd1 + 24);
    const auto *gpd1_27 = buffer.data(gpd1 + 27);
    const auto *gpd1_29 = buffer.data(gpd1 + 29);
    const auto *gpd1_40 = buffer.data(gpd1 + 40);
    const auto *gpd1_41 = buffer.data(gpd1 + 41);
    const auto *gpd1_48 = buffer.data(gpd1 + 48);
    const auto *gpd1_53 = buffer.data(gpd1 + 53);

    const auto *gpf_0 = buffer.data(gpf + 0);
    const auto *gpf_1 = buffer.data(gpf + 1);
    const auto *gpf_2 = buffer.data(gpf + 2);
    const auto *gpf_3 = buffer.data(gpf + 3);
    const auto *gpf_5 = buffer.data(gpf + 5);
    const auto *gpf_6 = buffer.data(gpf + 6);
    const auto *gpf_8 = buffer.data(gpf + 8);
    const auto *gpf_9 = buffer.data(gpf + 9);
    const auto *gpf_10 = buffer.data(gpf + 10);
    const auto *gpf_12 = buffer.data(gpf + 12);
    const auto *gpf_13 = buffer.data(gpf + 13);
    const auto *gpf_15 = buffer.data(gpf + 15);
    const auto *gpf_16 = buffer.data(gpf + 16);
    const auto *gpf_19 = buffer.data(gpf + 19);
    const auto *gpf_20 = buffer.data(gpf + 20);
    const auto *gpf_22 = buffer.data(gpf + 22);
    const auto *gpf_23 = buffer.data(gpf + 23);
    const auto *gpf_25 = buffer.data(gpf + 25);
    const auto *gpf_26 = buffer.data(gpf + 26);
    const auto *gpf_28 = buffer.data(gpf + 28);
    const auto *gpf_29 = buffer.data(gpf + 29);
    const auto *gpf_30 = buffer.data(gpf + 30);
    const auto *gpf_31 = buffer.data(gpf + 31);
    const auto *gpf_33 = buffer.data(gpf + 33);
    const auto *gpf_36 = buffer.data(gpf + 36);
    const auto *gpf_37 = buffer.data(gpf + 37);
    const auto *gpf_38 = buffer.data(gpf + 38);
    const auto *gpf_39 = buffer.data(gpf + 39);
    const auto *gpf_40 = buffer.data(gpf + 40);
    const auto *gpf_41 = buffer.data(gpf + 41);
    const auto *gpf_42 = buffer.data(gpf + 42);
    const auto *gpf_43 = buffer.data(gpf + 43);
    const auto *gpf_46 = buffer.data(gpf + 46);
    const auto *gpf_47 = buffer.data(gpf + 47);
    const auto *gpf_48 = buffer.data(gpf + 48);
    const auto *gpf_49 = buffer.data(gpf + 49);
    const auto *gpf_50 = buffer.data(gpf + 50);
    const auto *gpf_51 = buffer.data(gpf + 51);
    const auto *gpf_53 = buffer.data(gpf + 53);
    const auto *gpf_56 = buffer.data(gpf + 56);
    const auto *gpf_58 = buffer.data(gpf + 58);
    const auto *gpf_59 = buffer.data(gpf + 59);
    const auto *gpf_60 = buffer.data(gpf + 60);
    const auto *gpf_62 = buffer.data(gpf + 62);
    const auto *gpf_65 = buffer.data(gpf + 65);
    const auto *gpf_67 = buffer.data(gpf + 67);
    const auto *gpf_68 = buffer.data(gpf + 68);
    const auto *gpf_69 = buffer.data(gpf + 69);
    const auto *gpf_70 = buffer.data(gpf + 70);
    const auto *gpf_72 = buffer.data(gpf + 72);
    const auto *gpf_75 = buffer.data(gpf + 75);
    const auto *gpf_76 = buffer.data(gpf + 76);
    const auto *gpf_77 = buffer.data(gpf + 77);
    const auto *gpf_79 = buffer.data(gpf + 79);
    const auto *gpf_80 = buffer.data(gpf + 80);
    const auto *gpf_81 = buffer.data(gpf + 81);
    const auto *gpf_82 = buffer.data(gpf + 82);
    const auto *gpf_85 = buffer.data(gpf + 85);
    const auto *gpf_86 = buffer.data(gpf + 86);
    const auto *gpf_87 = buffer.data(gpf + 87);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, fpf_0, gsf_0, gpd0_0, \
                         gpd1_0, gpf_0, gpf_1, gpf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fpf_0[k]
                 + f_1 * gsf_0[k]
                 + f_2 * gpd0_0[k]
                 - f_3 * gpd1_0[k]
                 + f_4 * pc_x[k] * gpf_0[k];

        t_1[k] = f_4 * pc_y[k] * gpf_0[k];

        t_2[k] = f_4 * pc_z[k] * gpf_0[k];

        t_3[k] = f_5 * gpd0_0[k]
                 - f_6 * gpd1_0[k]
                 + f_4 * pc_y[k] * gpf_1[k];

        t_4[k] = f_4 * pc_y[k] * gpf_2[k];

        t_5[k] = f_5 * gpd0_0[k]
                 - f_6 * gpd1_0[k]
                 + f_4 * pc_z[k] * gpf_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pc_x, pc_y, pc_z, fpf_6, fpf_9, gsf_6, gsf_9, \
                         gpf_3, gpf_5, gpf_6, gpf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * fpf_6[k]
                 + f_1 * gsf_6[k]
                 + f_4 * pc_x[k] * gpf_6[k];

        t_7[k] = f_4 * pc_z[k] * gpf_3[k];

        t_8[k] = f_4 * pc_y[k] * gpf_5[k];

        t_9[k] = f_0 * fpf_9[k]
                 + f_1 * gsf_9[k]
                 + f_4 * pc_x[k] * gpf_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pc_y, pc_z, gpd0_3, gpd0_5, gpd1_3, \
                         gpd1_5, gpf_6, gpf_8, gpf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_2 * gpd0_3[k]
                  - f_3 * gpd1_3[k]
                  + f_4 * pc_y[k] * gpf_6[k];

        t_11[k] = f_4 * pc_z[k] * gpf_6[k];

        t_12[k] = f_5 * gpd0_5[k]
                  - f_6 * gpd1_5[k]
                  + f_4 * pc_y[k] * gpf_8[k];

        t_13[k] = f_4 * pc_y[k] * gpf_9[k];

        t_14[k] = f_2 * gpd0_5[k]
                  - f_3 * gpd1_5[k]
                  + f_4 * pc_z[k] * gpf_9[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pb_y, pc_y, pc_z, gsg0_0, gsg0_3, gsf_0, \
                         gsf_1, gsg1_0, gsg1_3, gpf_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pb_y[k] * gsg0_0[k]
                  - f_7 * pc_y[k] * gsg1_0[k];

        t_16[k] = f_1 * gsf_0[k]
                  + f_4 * pc_y[k] * gpf_10[k];

        t_17[k] = f_4 * pc_z[k] * gpf_10[k];

        t_18[k] = pb_y[k] * gsg0_3[k]
                  + f_8 * gsf_1[k]
                  - f_7 * pc_y[k] * gsg1_3[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pb_y, pc_x, pc_y, pc_z, fpf_16, gsg0_5, \
                         gsf_2, gsg1_5, gpf_12, gpf_13, gpf_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_1 * gsf_2[k]
                  + f_4 * pc_y[k] * gpf_12[k];

        t_20[k] = pb_y[k] * gsg0_5[k]
                  - f_7 * pc_y[k] * gsg1_5[k];

        t_21[k] = f_0 * fpf_16[k]
                  + f_4 * pc_x[k] * gpf_16[k];

        t_22[k] = f_4 * pc_z[k] * gpf_13[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pb_y, pc_x, pc_y, pc_z, fpf_19, gsg0_10, \
                         gsf_5, gsf_6, gsg1_10, gpf_15, gpf_16, \
                         gpf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * gsf_5[k]
                  + f_4 * pc_y[k] * gpf_15[k];

        t_24[k] = f_0 * fpf_19[k]
                  + f_4 * pc_x[k] * gpf_19[k];

        t_25[k] = pb_y[k] * gsg0_10[k]
                  + f_0 * gsf_6[k]
                  - f_7 * pc_y[k] * gsg1_10[k];

        t_26[k] = f_4 * pc_z[k] * gpf_16[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pb_y, pc_y, gsg0_12, gsg0_14, gsf_8, gsf_9, \
                         gsg1_12, gsg1_14, gpf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pb_y[k] * gsg0_12[k]
                  + f_8 * gsf_8[k]
                  - f_7 * pc_y[k] * gsg1_12[k];

        t_28[k] = f_1 * gsf_9[k]
                  + f_4 * pc_y[k] * gpf_19[k];

        t_29[k] = pb_y[k] * gsg0_14[k]
                  - f_7 * pc_y[k] * gsg1_14[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pb_z, pc_y, pc_z, gsg0_0, gsg0_3, \
                         gsf_0, gsg1_0, gsg1_3, gpf_20, gpf_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_z[k] * gsg0_0[k]
                  - f_7 * pc_z[k] * gsg1_0[k];

        t_31[k] = f_4 * pc_y[k] * gpf_20[k];

        t_32[k] = f_1 * gsf_0[k]
                  + f_4 * pc_z[k] * gpf_20[k];

        t_33[k] = pb_z[k] * gsg0_3[k]
                  - f_7 * pc_z[k] * gsg1_3[k];

        t_34[k] = f_4 * pc_y[k] * gpf_22[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pb_z, pc_x, pc_y, pc_z, fpf_26, gsg0_5, \
                         gsf_2, gsf_3, gsg1_5, gpf_23, gpf_25, gpf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pb_z[k] * gsg0_5[k]
                  + f_8 * gsf_2[k]
                  - f_7 * pc_z[k] * gsg1_5[k];

        t_36[k] = f_0 * fpf_26[k]
                  + f_4 * pc_x[k] * gpf_26[k];

        t_37[k] = f_1 * gsf_3[k]
                  + f_4 * pc_z[k] * gpf_23[k];

        t_38[k] = f_4 * pc_y[k] * gpf_25[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_z, pc_x, pc_z, fpf_29, gsg0_10, gsf_6, gsg1_10, \
                         gpf_26, gpf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_0 * fpf_29[k]
                  + f_4 * pc_x[k] * gpf_29[k];

        t_40[k] = pb_z[k] * gsg0_10[k]
                  - f_7 * pc_z[k] * gsg1_10[k];

        t_41[k] = f_1 * gsf_6[k]
                  + f_4 * pc_z[k] * gpf_26[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pb_z, pc_y, pc_z, gsg0_14, gsf_9, gsg1_14, gpd0_17, \
                         gpd1_17, gpf_28, gpf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_5 * gpd0_17[k]
                  - f_6 * gpd1_17[k]
                  + f_4 * pc_y[k] * gpf_28[k];

        t_43[k] = f_4 * pc_y[k] * gpf_29[k];

        t_44[k] = pb_z[k] * gsg0_14[k]
                  + f_0 * gsf_9[k]
                  - f_7 * pc_z[k] * gsg1_14[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, pa_y, pc_y, pc_z, fpg0_0, fpg0_3, \
                         fpf_0, fpf_1, fpg1_0, fpg1_3, gpf_30, gpf_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = pa_y[k] * fpg0_0[k]
                  - f_7 * pc_y[k] * fpg1_0[k];

        t_46[k] = f_1 * fpf_0[k]
                  + f_4 * pc_y[k] * gpf_30[k];

        t_47[k] = f_4 * pc_z[k] * gpf_30[k];

        t_48[k] = pa_y[k] * fpg0_3[k]
                  + f_8 * fpf_1[k]
                  - f_7 * pc_y[k] * fpg1_3[k];

        t_49[k] = f_4 * pc_z[k] * gpf_31[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pa_y, pc_x, pc_y, pc_z, fpg0_5, fpf_36, fpg1_5, \
                         gsf_16, gpf_33, gpf_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_y[k] * fpg0_5[k]
                  - f_7 * pc_y[k] * fpg1_5[k];

        t_51[k] = f_9 * fpf_36[k]
                  + f_1 * gsf_16[k]
                  + f_4 * pc_x[k] * gpf_36[k];

        t_52[k] = f_4 * pc_z[k] * gpf_33[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pa_y, pc_x, pc_y, fpg0_9, fpf_6, fpf_38, fpg1_9, \
                         gsf_18, gpd0_21, gpd1_21, gpf_36, gpf_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_9 * fpf_38[k]
                  + f_1 * gsf_18[k]
                  + f_4 * pc_x[k] * gpf_38[k];

        t_54[k] = pa_y[k] * fpg0_9[k]
                  - f_7 * pc_y[k] * fpg1_9[k];

        t_55[k] = f_1 * fpf_6[k]
                  + f_2 * gpd0_21[k]
                  - f_3 * gpd1_21[k]
                  + f_4 * pc_y[k] * gpf_36[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_y, pc_y, pc_z, fpg0_14, fpf_9, fpg1_14, \
                         gpd0_21, gpd1_21, gpf_36, gpf_37, gpf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_4 * pc_z[k] * gpf_36[k];

        t_57[k] = f_5 * gpd0_21[k]
                  - f_6 * gpd1_21[k]
                  + f_4 * pc_z[k] * gpf_37[k];

        t_58[k] = f_1 * fpf_9[k]
                  + f_4 * pc_y[k] * gpf_39[k];

        t_59[k] = pa_y[k] * fpg0_14[k]
                  - f_7 * pc_y[k] * fpg1_14[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, pc_x, pc_y, pc_z, fpf_10, fpf_40, gsf_10, gpd0_24, \
                         gpd1_24, gpf_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_9 * fpf_40[k]
                  + f_2 * gpd0_24[k]
                  - f_3 * gpd1_24[k]
                  + f_4 * pc_x[k] * gpf_40[k];

        t_61[k] = f_1 * fpf_10[k]
                  + f_1 * gsf_10[k]
                  + f_4 * pc_y[k] * gpf_40[k];

        t_62[k] = f_4 * pc_z[k] * gpf_40[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pc_x, pc_z, fpf_43, fpf_46, gpd0_24, gpd0_27, \
                         gpd1_24, gpd1_27, gpf_41, gpf_42, gpf_43, \
                         gpf_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_9 * fpf_43[k]
                  + f_5 * gpd0_27[k]
                  - f_6 * gpd1_27[k]
                  + f_4 * pc_x[k] * gpf_43[k];

        t_64[k] = f_4 * pc_z[k] * gpf_41[k];

        t_65[k] = f_5 * gpd0_24[k]
                  - f_6 * gpd1_24[k]
                  + f_4 * pc_z[k] * gpf_42[k];

        t_66[k] = f_9 * fpf_46[k]
                  + f_4 * pc_x[k] * gpf_46[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_x, pc_x, pc_z, dpg0_70, dpg1_70, fpg0_70, \
                         fpf_48, fpf_49, fpg1_70, gpf_43, gpf_48, \
                         gpf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_4 * pc_z[k] * gpf_43[k];

        t_68[k] = f_9 * fpf_48[k]
                  + f_4 * pc_x[k] * gpf_48[k];

        t_69[k] = f_9 * fpf_49[k]
                  + f_4 * pc_x[k] * gpf_49[k];

        t_70[k] = f_10 * dpg0_70[k]
                  - f_11 * dpg1_70[k]
                  + pa_x[k] * fpg0_70[k]
                  - f_7 * pc_x[k] * fpg1_70[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pc_y, pc_z, fpf_19, gsf_19, gpd0_27, gpd0_29, \
                         gpd1_27, gpd1_29, gpf_46, gpf_47, gpf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_4 * pc_z[k] * gpf_46[k];

        t_72[k] = f_5 * gpd0_27[k]
                  - f_6 * gpd1_27[k]
                  + f_4 * pc_z[k] * gpf_47[k];

        t_73[k] = f_1 * fpf_19[k]
                  + f_1 * gsf_19[k]
                  + f_4 * pc_y[k] * gpf_49[k];

        t_74[k] = f_2 * gpd0_29[k]
                  - f_3 * gpd1_29[k]
                  + f_4 * pc_z[k] * gpf_49[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pa_y, pb_z, pc_y, pc_z, fpg0_30, fpf_20, \
                         fpg1_30, gsg0_18, gsf_10, gsg1_18, gpf_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = pa_y[k] * fpg0_30[k]
                  - f_7 * pc_y[k] * fpg1_30[k];

        t_76[k] = f_1 * fpf_20[k]
                  + f_4 * pc_y[k] * gpf_50[k];

        t_77[k] = f_1 * gsf_10[k]
                  + f_4 * pc_z[k] * gpf_50[k];

        t_78[k] = pb_z[k] * gsg0_18[k]
                  - f_7 * pc_z[k] * gsg1_18[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pa_y, pc_x, pc_y, pc_z, fpg0_35, fpf_56, \
                         fpg1_35, gsf_11, gsf_13, gpf_51, gpf_53, \
                         gpf_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_1 * gsf_11[k]
                  + f_4 * pc_z[k] * gpf_51[k];

        t_80[k] = pa_y[k] * fpg0_35[k]
                  - f_7 * pc_y[k] * fpg1_35[k];

        t_81[k] = f_9 * fpf_56[k]
                  + f_4 * pc_x[k] * gpf_56[k];

        t_82[k] = f_1 * gsf_13[k]
                  + f_4 * pc_z[k] * gpf_53[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pb_z, pc_x, pc_z, fpf_58, fpf_59, gsg0_25, \
                         gsf_16, gsg1_25, gpf_56, gpf_58, gpf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_9 * fpf_58[k]
                  + f_4 * pc_x[k] * gpf_58[k];

        t_84[k] = f_9 * fpf_59[k]
                  + f_4 * pc_x[k] * gpf_59[k];

        t_85[k] = pb_z[k] * gsg0_25[k]
                  - f_7 * pc_z[k] * gsg1_25[k];

        t_86[k] = f_1 * gsf_16[k]
                  + f_4 * pc_z[k] * gpf_56[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pa_y, pb_z, pc_y, pc_z, fpg0_44, fpf_29, fpg1_44, \
                         gsg0_27, gsf_17, gsg1_27, gpf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pb_z[k] * gsg0_27[k]
                  + f_8 * gsf_17[k]
                  - f_7 * pc_z[k] * gsg1_27[k];

        t_88[k] = f_1 * fpf_29[k]
                  + f_4 * pc_y[k] * gpf_59[k];

        t_89[k] = pa_y[k] * fpg0_44[k]
                  - f_7 * pc_y[k] * fpg1_44[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, pa_z, pc_y, pc_z, fpg0_0, fpg0_3, \
                         fpf_0, fpg1_0, fpg1_3, gpf_60, gpf_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = pa_z[k] * fpg0_0[k]
                  - f_7 * pc_z[k] * fpg1_0[k];

        t_91[k] = f_4 * pc_y[k] * gpf_60[k];

        t_92[k] = f_1 * fpf_0[k]
                  + f_4 * pc_z[k] * gpf_60[k];

        t_93[k] = pa_z[k] * fpg0_3[k]
                  - f_7 * pc_z[k] * fpg1_3[k];

        t_94[k] = f_4 * pc_y[k] * gpf_62[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, pa_z, pc_x, pc_z, fpg0_5, fpg0_6, fpf_2, fpf_67, \
                         fpg1_5, fpg1_6, gsf_27, gpf_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = pa_z[k] * fpg0_5[k]
                  + f_8 * fpf_2[k]
                  - f_7 * pc_z[k] * fpg1_5[k];

        t_96[k] = pa_z[k] * fpg0_6[k]
                  - f_7 * pc_z[k] * fpg1_6[k];

        t_97[k] = f_9 * fpf_67[k]
                  + f_1 * gsf_27[k]
                  + f_4 * pc_x[k] * gpf_67[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, pa_z, pc_x, pc_y, pc_z, fpg0_10, fpf_69, fpg1_10, \
                         gsf_29, gpf_65, gpf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_4 * pc_y[k] * gpf_65[k];

        t_99[k] = f_9 * fpf_69[k]
                  + f_1 * gsf_29[k]
                  + f_4 * pc_x[k] * gpf_69[k];

        t_100[k] = pa_z[k] * fpg0_10[k]
                   - f_7 * pc_z[k] * fpg1_10[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pc_y, pc_z, fpf_9, gpd0_40, gpd0_41, \
                         gpd1_40, gpd1_41, gpf_67, gpf_68, gpf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_12 * gpd0_40[k]
                   - f_13 * gpd1_40[k]
                   + f_4 * pc_y[k] * gpf_67[k];

        t_102[k] = f_5 * gpd0_41[k]
                   - f_6 * gpd1_41[k]
                   + f_4 * pc_y[k] * gpf_68[k];

        t_103[k] = f_4 * pc_y[k] * gpf_69[k];

        t_104[k] = f_1 * fpf_9[k]
                   + f_2 * gpd0_41[k]
                   - f_3 * gpd1_41[k]
                   + f_4 * pc_z[k] * gpf_69[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pa_z, pc_y, pc_z, fpg0_15, fpg0_18, \
                         fpf_10, fpg1_15, fpg1_18, gsf_20, gpf_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = pa_z[k] * fpg0_15[k]
                   - f_7 * pc_z[k] * fpg1_15[k];

        t_106[k] = f_1 * gsf_20[k]
                   + f_4 * pc_y[k] * gpf_70[k];

        t_107[k] = f_1 * fpf_10[k]
                   + f_4 * pc_z[k] * gpf_70[k];

        t_108[k] = pa_z[k] * fpg0_18[k]
                   - f_7 * pc_z[k] * fpg1_18[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pb_y, pc_x, pc_y, fpf_76, fpf_77, \
                         gsg0_35, gsf_22, gsg1_35, gpf_72, gpf_76, \
                         gpf_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_1 * gsf_22[k]
                   + f_4 * pc_y[k] * gpf_72[k];

        t_110[k] = pb_y[k] * gsg0_35[k]
                   - f_7 * pc_y[k] * gsg1_35[k];

        t_111[k] = f_9 * fpf_76[k]
                   + f_4 * pc_x[k] * gpf_76[k];

        t_112[k] = f_9 * fpf_77[k]
                   + f_4 * pc_x[k] * gpf_77[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, pa_z, pc_x, pc_y, pc_z, fpg0_25, fpf_79, \
                         fpg1_25, gsf_25, gpf_75, gpf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_1 * gsf_25[k]
                   + f_4 * pc_y[k] * gpf_75[k];

        t_114[k] = f_9 * fpf_79[k]
                   + f_4 * pc_x[k] * gpf_79[k];

        t_115[k] = pa_z[k] * fpg0_25[k]
                   - f_7 * pc_z[k] * fpg1_25[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pb_y, pc_y, gsg0_41, gsg0_42, gsg0_44, \
                         gsf_27, gsf_28, gsf_29, gsg1_41, gsg1_42, gsg1_44, \
                         gpf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = pb_y[k] * gsg0_41[k]
                   + f_9 * gsf_27[k]
                   - f_7 * pc_y[k] * gsg1_41[k];

        t_117[k] = pb_y[k] * gsg0_42[k]
                   + f_8 * gsf_28[k]
                   - f_7 * pc_y[k] * gsg1_42[k];

        t_118[k] = f_1 * gsf_29[k]
                   + f_4 * pc_y[k] * gpf_79[k];

        t_119[k] = pb_y[k] * gsg0_44[k]
                   - f_7 * pc_y[k] * gsg1_44[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, t_124, pc_x, pc_y, pc_z, fpf_20, fpf_80, \
                         gsf_20, gpd0_48, gpd1_48, gpf_80, gpf_81, \
                         gpf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_9 * fpf_80[k]
                   + f_2 * gpd0_48[k]
                   - f_3 * gpd1_48[k]
                   + f_4 * pc_x[k] * gpf_80[k];

        t_121[k] = f_4 * pc_y[k] * gpf_80[k];

        t_122[k] = f_1 * fpf_20[k]
                   + f_1 * gsf_20[k]
                   + f_4 * pc_z[k] * gpf_80[k];

        t_123[k] = f_5 * gpd0_48[k]
                   - f_6 * gpd1_48[k]
                   + f_4 * pc_y[k] * gpf_81[k];

        t_124[k] = f_4 * pc_y[k] * gpf_82[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pc_x, pc_y, fpf_85, fpf_86, fpf_87, \
                         gpd0_53, gpd1_53, gpf_85, gpf_86, gpf_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_9 * fpf_85[k]
                   + f_5 * gpd0_53[k]
                   - f_6 * gpd1_53[k]
                   + f_4 * pc_x[k] * gpf_85[k];

        t_126[k] = f_9 * fpf_86[k]
                   + f_4 * pc_x[k] * gpf_86[k];

        t_127[k] = f_9 * fpf_87[k]
                   + f_4 * pc_x[k] * gpf_87[k];

        t_128[k] = f_4 * pc_y[k] * gpf_85[k];
    }
}

static auto
compute_prim_gpg_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t dpg0, const size_t dpg1,
                                                          const size_t fpg0, const size_t fpf,
                                                          const size_t fpg1, const size_t gsg0,
                                                          const size_t gsf, const size_t gsg1,
                                                          const size_t gpd0, const size_t gpd1,
                                                          const size_t gpf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 1.5 / gamma;
    const auto f_3 = 1.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = gamma / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 1.5 / q;
    const auto f_10 = 1.0 / p;
    const auto f_11 = gamma / (p * q);
    const auto f_12 = 1.0 / gamma;
    const auto f_13 = p / (gamma * q);
    const auto f_14 = 0.5 / p;
    const auto f_15 = 0.5 * gamma / (p * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dpg0_0 = buffer.data(dpg0 + 0);
    const auto *dpg0_134 = buffer.data(dpg0 + 134);
    const auto *dpg0_160 = buffer.data(dpg0 + 160);
    const auto *dpg0_207 = buffer.data(dpg0 + 207);

    const auto *dpg1_0 = buffer.data(dpg1 + 0);
    const auto *dpg1_134 = buffer.data(dpg1 + 134);
    const auto *dpg1_160 = buffer.data(dpg1 + 160);
    const auto *dpg1_207 = buffer.data(dpg1 + 207);

    const auto *fpg0_45 = buffer.data(fpg0 + 45);
    const auto *fpg0_48 = buffer.data(fpg0 + 48);
    const auto *fpg0_51 = buffer.data(fpg0 + 51);
    const auto *fpg0_55 = buffer.data(fpg0 + 55);
    const auto *fpg0_61 = buffer.data(fpg0 + 61);
    const auto *fpg0_63 = buffer.data(fpg0 + 63);
    const auto *fpg0_70 = buffer.data(fpg0 + 70);
    const auto *fpg0_90 = buffer.data(fpg0 + 90);
    const auto *fpg0_95 = buffer.data(fpg0 + 95);
    const auto *fpg0_99 = buffer.data(fpg0 + 99);
    const auto *fpg0_104 = buffer.data(fpg0 + 104);
    const auto *fpg0_120 = buffer.data(fpg0 + 120);
    const auto *fpg0_122 = buffer.data(fpg0 + 122);
    const auto *fpg0_125 = buffer.data(fpg0 + 125);
    const auto *fpg0_134 = buffer.data(fpg0 + 134);
    const auto *fpg0_160 = buffer.data(fpg0 + 160);
    const auto *fpg0_207 = buffer.data(fpg0 + 207);

    const auto *fpf_30 = buffer.data(fpf + 30);
    const auto *fpf_36 = buffer.data(fpf + 36);
    const auto *fpf_39 = buffer.data(fpf + 39);
    const auto *fpf_40 = buffer.data(fpf + 40);
    const auto *fpf_46 = buffer.data(fpf + 46);
    const auto *fpf_49 = buffer.data(fpf + 49);
    const auto *fpf_50 = buffer.data(fpf + 50);
    const auto *fpf_56 = buffer.data(fpf + 56);
    const auto *fpf_59 = buffer.data(fpf + 59);
    const auto *fpf_60 = buffer.data(fpf + 60);
    const auto *fpf_62 = buffer.data(fpf + 62);
    const auto *fpf_68 = buffer.data(fpf + 68);
    const auto *fpf_69 = buffer.data(fpf + 69);
    const auto *fpf_70 = buffer.data(fpf + 70);
    const auto *fpf_72 = buffer.data(fpf + 72);
    const auto *fpf_79 = buffer.data(fpf + 79);
    const auto *fpf_80 = buffer.data(fpf + 80);
    const auto *fpf_82 = buffer.data(fpf + 82);
    const auto *fpf_86 = buffer.data(fpf + 86);
    const auto *fpf_88 = buffer.data(fpf + 88);
    const auto *fpf_89 = buffer.data(fpf + 89);
    const auto *fpf_93 = buffer.data(fpf + 93);
    const auto *fpf_96 = buffer.data(fpf + 96);
    const auto *fpf_98 = buffer.data(fpf + 98);
    const auto *fpf_99 = buffer.data(fpf + 99);
    const auto *fpf_100 = buffer.data(fpf + 100);
    const auto *fpf_103 = buffer.data(fpf + 103);
    const auto *fpf_106 = buffer.data(fpf + 106);
    const auto *fpf_108 = buffer.data(fpf + 108);
    const auto *fpf_109 = buffer.data(fpf + 109);
    const auto *fpf_116 = buffer.data(fpf + 116);
    const auto *fpf_118 = buffer.data(fpf + 118);
    const auto *fpf_119 = buffer.data(fpf + 119);
    const auto *fpf_127 = buffer.data(fpf + 127);
    const auto *fpf_128 = buffer.data(fpf + 128);
    const auto *fpf_130 = buffer.data(fpf + 130);
    const auto *fpf_135 = buffer.data(fpf + 135);
    const auto *fpf_136 = buffer.data(fpf + 136);
    const auto *fpf_137 = buffer.data(fpf + 137);
    const auto *fpf_138 = buffer.data(fpf + 138);
    const auto *fpf_139 = buffer.data(fpf + 139);
    const auto *fpf_143 = buffer.data(fpf + 143);
    const auto *fpf_146 = buffer.data(fpf + 146);
    const auto *fpf_147 = buffer.data(fpf + 147);
    const auto *fpf_148 = buffer.data(fpf + 148);
    const auto *fpf_149 = buffer.data(fpf + 149);
    const auto *fpf_155 = buffer.data(fpf + 155);
    const auto *fpf_156 = buffer.data(fpf + 156);
    const auto *fpf_157 = buffer.data(fpf + 157);
    const auto *fpf_159 = buffer.data(fpf + 159);

    const auto *fpg1_45 = buffer.data(fpg1 + 45);
    const auto *fpg1_48 = buffer.data(fpg1 + 48);
    const auto *fpg1_51 = buffer.data(fpg1 + 51);
    const auto *fpg1_55 = buffer.data(fpg1 + 55);
    const auto *fpg1_61 = buffer.data(fpg1 + 61);
    const auto *fpg1_63 = buffer.data(fpg1 + 63);
    const auto *fpg1_70 = buffer.data(fpg1 + 70);
    const auto *fpg1_90 = buffer.data(fpg1 + 90);
    const auto *fpg1_95 = buffer.data(fpg1 + 95);
    const auto *fpg1_99 = buffer.data(fpg1 + 99);
    const auto *fpg1_104 = buffer.data(fpg1 + 104);
    const auto *fpg1_120 = buffer.data(fpg1 + 120);
    const auto *fpg1_122 = buffer.data(fpg1 + 122);
    const auto *fpg1_125 = buffer.data(fpg1 + 125);
    const auto *fpg1_134 = buffer.data(fpg1 + 134);
    const auto *fpg1_160 = buffer.data(fpg1 + 160);
    const auto *fpg1_207 = buffer.data(fpg1 + 207);

    const auto *gsg0_45 = buffer.data(gsg0 + 45);
    const auto *gsg0_48 = buffer.data(gsg0 + 48);
    const auto *gsg0_50 = buffer.data(gsg0 + 50);
    const auto *gsg0_55 = buffer.data(gsg0 + 55);
    const auto *gsg0_57 = buffer.data(gsg0 + 57);
    const auto *gsg0_59 = buffer.data(gsg0 + 59);
    const auto *gsg0_75 = buffer.data(gsg0 + 75);
    const auto *gsg0_78 = buffer.data(gsg0 + 78);

    const auto *gsf_30 = buffer.data(gsf + 30);
    const auto *gsf_31 = buffer.data(gsf + 31);
    const auto *gsf_32 = buffer.data(gsf + 32);
    const auto *gsf_33 = buffer.data(gsf + 33);
    const auto *gsf_36 = buffer.data(gsf + 36);
    const auto *gsf_37 = buffer.data(gsf + 37);
    const auto *gsf_38 = buffer.data(gsf + 38);
    const auto *gsf_39 = buffer.data(gsf + 39);
    const auto *gsf_42 = buffer.data(gsf + 42);
    const auto *gsf_46 = buffer.data(gsf + 46);
    const auto *gsf_47 = buffer.data(gsf + 47);
    const auto *gsf_48 = buffer.data(gsf + 48);
    const auto *gsf_49 = buffer.data(gsf + 49);
    const auto *gsf_50 = buffer.data(gsf + 50);
    const auto *gsf_51 = buffer.data(gsf + 51);
    const auto *gsf_55 = buffer.data(gsf + 55);
    const auto *gsf_56 = buffer.data(gsf + 56);
    const auto *gsf_57 = buffer.data(gsf + 57);
    const auto *gsf_59 = buffer.data(gsf + 59);

    const auto *gsg1_45 = buffer.data(gsg1 + 45);
    const auto *gsg1_48 = buffer.data(gsg1 + 48);
    const auto *gsg1_50 = buffer.data(gsg1 + 50);
    const auto *gsg1_55 = buffer.data(gsg1 + 55);
    const auto *gsg1_57 = buffer.data(gsg1 + 57);
    const auto *gsg1_59 = buffer.data(gsg1 + 59);
    const auto *gsg1_75 = buffer.data(gsg1 + 75);
    const auto *gsg1_78 = buffer.data(gsg1 + 78);

    const auto *gpd0_51 = buffer.data(gpd0 + 51);
    const auto *gpd0_52 = buffer.data(gpd0 + 52);
    const auto *gpd0_53 = buffer.data(gpd0 + 53);
    const auto *gpd0_54 = buffer.data(gpd0 + 54);
    const auto *gpd0_57 = buffer.data(gpd0 + 57);
    const auto *gpd0_59 = buffer.data(gpd0 + 59);
    const auto *gpd0_60 = buffer.data(gpd0 + 60);
    const auto *gpd0_63 = buffer.data(gpd0 + 63);
    const auto *gpd0_65 = buffer.data(gpd0 + 65);
    const auto *gpd0_77 = buffer.data(gpd0 + 77);
    const auto *gpd0_78 = buffer.data(gpd0 + 78);
    const auto *gpd0_83 = buffer.data(gpd0 + 83);
    const auto *gpd0_87 = buffer.data(gpd0 + 87);
    const auto *gpd0_89 = buffer.data(gpd0 + 89);
    const auto *gpd0_90 = buffer.data(gpd0 + 90);
    const auto *gpd0_93 = buffer.data(gpd0 + 93);
    const auto *gpd0_94 = buffer.data(gpd0 + 94);
    const auto *gpd0_95 = buffer.data(gpd0 + 95);

    const auto *gpd1_51 = buffer.data(gpd1 + 51);
    const auto *gpd1_52 = buffer.data(gpd1 + 52);
    const auto *gpd1_53 = buffer.data(gpd1 + 53);
    const auto *gpd1_54 = buffer.data(gpd1 + 54);
    const auto *gpd1_57 = buffer.data(gpd1 + 57);
    const auto *gpd1_59 = buffer.data(gpd1 + 59);
    const auto *gpd1_60 = buffer.data(gpd1 + 60);
    const auto *gpd1_63 = buffer.data(gpd1 + 63);
    const auto *gpd1_65 = buffer.data(gpd1 + 65);
    const auto *gpd1_77 = buffer.data(gpd1 + 77);
    const auto *gpd1_78 = buffer.data(gpd1 + 78);
    const auto *gpd1_83 = buffer.data(gpd1 + 83);
    const auto *gpd1_87 = buffer.data(gpd1 + 87);
    const auto *gpd1_89 = buffer.data(gpd1 + 89);
    const auto *gpd1_90 = buffer.data(gpd1 + 90);
    const auto *gpd1_93 = buffer.data(gpd1 + 93);
    const auto *gpd1_94 = buffer.data(gpd1 + 94);
    const auto *gpd1_95 = buffer.data(gpd1 + 95);

    const auto *gpf_86 = buffer.data(gpf + 86);
    const auto *gpf_87 = buffer.data(gpf + 87);
    const auto *gpf_88 = buffer.data(gpf + 88);
    const auto *gpf_89 = buffer.data(gpf + 89);
    const auto *gpf_90 = buffer.data(gpf + 90);
    const auto *gpf_91 = buffer.data(gpf + 91);
    const auto *gpf_92 = buffer.data(gpf + 92);
    const auto *gpf_93 = buffer.data(gpf + 93);
    const auto *gpf_96 = buffer.data(gpf + 96);
    const auto *gpf_97 = buffer.data(gpf + 97);
    const auto *gpf_98 = buffer.data(gpf + 98);
    const auto *gpf_99 = buffer.data(gpf + 99);
    const auto *gpf_100 = buffer.data(gpf + 100);
    const auto *gpf_101 = buffer.data(gpf + 101);
    const auto *gpf_102 = buffer.data(gpf + 102);
    const auto *gpf_103 = buffer.data(gpf + 103);
    const auto *gpf_106 = buffer.data(gpf + 106);
    const auto *gpf_107 = buffer.data(gpf + 107);
    const auto *gpf_108 = buffer.data(gpf + 108);
    const auto *gpf_109 = buffer.data(gpf + 109);
    const auto *gpf_110 = buffer.data(gpf + 110);
    const auto *gpf_111 = buffer.data(gpf + 111);
    const auto *gpf_113 = buffer.data(gpf + 113);
    const auto *gpf_116 = buffer.data(gpf + 116);
    const auto *gpf_118 = buffer.data(gpf + 118);
    const auto *gpf_119 = buffer.data(gpf + 119);
    const auto *gpf_120 = buffer.data(gpf + 120);
    const auto *gpf_122 = buffer.data(gpf + 122);
    const auto *gpf_126 = buffer.data(gpf + 126);
    const auto *gpf_127 = buffer.data(gpf + 127);
    const auto *gpf_128 = buffer.data(gpf + 128);
    const auto *gpf_129 = buffer.data(gpf + 129);
    const auto *gpf_130 = buffer.data(gpf + 130);
    const auto *gpf_132 = buffer.data(gpf + 132);
    const auto *gpf_135 = buffer.data(gpf + 135);
    const auto *gpf_136 = buffer.data(gpf + 136);
    const auto *gpf_137 = buffer.data(gpf + 137);
    const auto *gpf_138 = buffer.data(gpf + 138);
    const auto *gpf_139 = buffer.data(gpf + 139);
    const auto *gpf_140 = buffer.data(gpf + 140);
    const auto *gpf_142 = buffer.data(gpf + 142);
    const auto *gpf_143 = buffer.data(gpf + 143);
    const auto *gpf_146 = buffer.data(gpf + 146);
    const auto *gpf_147 = buffer.data(gpf + 147);
    const auto *gpf_148 = buffer.data(gpf + 148);
    const auto *gpf_149 = buffer.data(gpf + 149);
    const auto *gpf_150 = buffer.data(gpf + 150);
    const auto *gpf_151 = buffer.data(gpf + 151);
    const auto *gpf_152 = buffer.data(gpf + 152);
    const auto *gpf_155 = buffer.data(gpf + 155);
    const auto *gpf_156 = buffer.data(gpf + 156);
    const auto *gpf_157 = buffer.data(gpf + 157);
    const auto *gpf_158 = buffer.data(gpf + 158);
    const auto *gpf_159 = buffer.data(gpf + 159);
    const auto *gpf_160 = buffer.data(gpf + 160);

#pragma omp simd aligned(t_129, t_130, t_131, pc_x, pc_y, fpf_89, gpd0_51, gpd0_52, gpd1_51, \
                         gpd1_52, gpf_86, gpf_87, gpf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_9 * fpf_89[k]
                   + f_4 * pc_x[k] * gpf_89[k];

        t_130[k] = f_2 * gpd0_51[k]
                   - f_3 * gpd1_51[k]
                   + f_4 * pc_y[k] * gpf_86[k];

        t_131[k] = f_12 * gpd0_52[k]
                   - f_13 * gpd1_52[k]
                   + f_4 * pc_y[k] * gpf_87[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, pa_x, pc_x, pc_y, dpg0_134, dpg1_134, fpg0_134, \
                         fpg1_134, gpd0_53, gpd1_53, gpf_88, gpf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_5 * gpd0_53[k]
                   - f_6 * gpd1_53[k]
                   + f_4 * pc_y[k] * gpf_88[k];

        t_133[k] = f_4 * pc_y[k] * gpf_89[k];

        t_134[k] = f_10 * dpg0_134[k]
                   - f_11 * dpg1_134[k]
                   + pa_x[k] * fpg0_134[k]
                   - f_7 * pc_x[k] * fpg1_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, pa_y, pc_y, pc_z, dpg0_0, dpg1_0, fpg0_45, \
                         fpf_30, fpg1_45, gpf_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_14 * dpg0_0[k]
                   - f_15 * dpg1_0[k]
                   + pa_y[k] * fpg0_45[k]
                   - f_7 * pc_y[k] * fpg1_45[k];

        t_136[k] = f_8 * fpf_30[k]
                   + f_4 * pc_y[k] * gpf_90[k];

        t_137[k] = f_4 * pc_z[k] * gpf_90[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pc_x, pc_z, fpf_93, gsf_33, gpd0_54, gpd0_57, \
                         gpd1_54, gpd1_57, gpf_91, gpf_92, gpf_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_8 * fpf_93[k]
                   + f_1 * gsf_33[k]
                   + f_5 * gpd0_57[k]
                   - f_6 * gpd1_57[k]
                   + f_4 * pc_x[k] * gpf_93[k];

        t_139[k] = f_4 * pc_z[k] * gpf_91[k];

        t_140[k] = f_5 * gpd0_54[k]
                   - f_6 * gpd1_54[k]
                   + f_4 * pc_z[k] * gpf_92[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, t_144, pc_x, pc_z, fpf_96, fpf_98, fpf_99, \
                         gsf_36, gsf_38, gsf_39, gpf_93, gpf_96, gpf_98, \
                         gpf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_8 * fpf_96[k]
                   + f_1 * gsf_36[k]
                   + f_4 * pc_x[k] * gpf_96[k];

        t_142[k] = f_4 * pc_z[k] * gpf_93[k];

        t_143[k] = f_8 * fpf_98[k]
                   + f_1 * gsf_38[k]
                   + f_4 * pc_x[k] * gpf_98[k];

        t_144[k] = f_8 * fpf_99[k]
                   + f_1 * gsf_39[k]
                   + f_4 * pc_x[k] * gpf_99[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, t_149, pc_y, pc_z, fpf_36, fpf_39, \
                         gpd0_57, gpd0_59, gpd1_57, gpd1_59, gpf_96, gpf_97, \
                         gpf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_8 * fpf_36[k]
                   + f_2 * gpd0_57[k]
                   - f_3 * gpd1_57[k]
                   + f_4 * pc_y[k] * gpf_96[k];

        t_146[k] = f_4 * pc_z[k] * gpf_96[k];

        t_147[k] = f_5 * gpd0_57[k]
                   - f_6 * gpd1_57[k]
                   + f_4 * pc_z[k] * gpf_97[k];

        t_148[k] = f_8 * fpf_39[k]
                   + f_4 * pc_y[k] * gpf_99[k];

        t_149[k] = f_2 * gpd0_59[k]
                   - f_3 * gpd1_59[k]
                   + f_4 * pc_z[k] * gpf_99[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pc_x, pc_y, pc_z, fpf_40, fpf_100, gsf_30, \
                         gpd0_60, gpd1_60, gpf_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_8 * fpf_100[k]
                   + f_2 * gpd0_60[k]
                   - f_3 * gpd1_60[k]
                   + f_4 * pc_x[k] * gpf_100[k];

        t_151[k] = f_8 * fpf_40[k]
                   + f_1 * gsf_30[k]
                   + f_4 * pc_y[k] * gpf_100[k];

        t_152[k] = f_4 * pc_z[k] * gpf_100[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, pc_x, pc_z, fpf_103, fpf_106, gpd0_60, \
                         gpd0_63, gpd1_60, gpd1_63, gpf_101, gpf_102, gpf_103, \
                         gpf_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_8 * fpf_103[k]
                   + f_5 * gpd0_63[k]
                   - f_6 * gpd1_63[k]
                   + f_4 * pc_x[k] * gpf_103[k];

        t_154[k] = f_4 * pc_z[k] * gpf_101[k];

        t_155[k] = f_5 * gpd0_60[k]
                   - f_6 * gpd1_60[k]
                   + f_4 * pc_z[k] * gpf_102[k];

        t_156[k] = f_8 * fpf_106[k]
                   + f_4 * pc_x[k] * gpf_106[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, pa_x, pc_x, pc_z, dpg0_160, dpg1_160, \
                         fpg0_160, fpf_108, fpf_109, fpg1_160, gpf_103, gpf_108, \
                         gpf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_4 * pc_z[k] * gpf_103[k];

        t_158[k] = f_8 * fpf_108[k]
                   + f_4 * pc_x[k] * gpf_108[k];

        t_159[k] = f_8 * fpf_109[k]
                   + f_4 * pc_x[k] * gpf_109[k];

        t_160[k] = f_14 * dpg0_160[k]
                   - f_15 * dpg1_160[k]
                   + pa_x[k] * fpg0_160[k]
                   - f_7 * pc_x[k] * fpg1_160[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, pc_y, pc_z, fpf_49, gsf_39, gpd0_63, \
                         gpd0_65, gpd1_63, gpd1_65, gpf_106, gpf_107, \
                         gpf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_4 * pc_z[k] * gpf_106[k];

        t_162[k] = f_5 * gpd0_63[k]
                   - f_6 * gpd1_63[k]
                   + f_4 * pc_z[k] * gpf_107[k];

        t_163[k] = f_8 * fpf_49[k]
                   + f_1 * gsf_39[k]
                   + f_4 * pc_y[k] * gpf_109[k];

        t_164[k] = f_2 * gpd0_65[k]
                   - f_3 * gpd1_65[k]
                   + f_4 * pc_z[k] * gpf_109[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, t_168, pb_z, pc_y, pc_z, fpf_50, gsg0_45, \
                         gsg0_48, gsf_30, gsg1_45, gsg1_48, gpf_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = pb_z[k] * gsg0_45[k]
                   - f_7 * pc_z[k] * gsg1_45[k];

        t_166[k] = f_8 * fpf_50[k]
                   + f_4 * pc_y[k] * gpf_110[k];

        t_167[k] = f_1 * gsf_30[k]
                   + f_4 * pc_z[k] * gpf_110[k];

        t_168[k] = pb_z[k] * gsg0_48[k]
                   - f_7 * pc_z[k] * gsg1_48[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, pb_z, pc_x, pc_z, fpf_116, gsg0_50, \
                         gsf_31, gsf_32, gsf_33, gsg1_50, gpf_111, gpf_113, \
                         gpf_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = f_1 * gsf_31[k]
                   + f_4 * pc_z[k] * gpf_111[k];

        t_170[k] = pb_z[k] * gsg0_50[k]
                   + f_8 * gsf_32[k]
                   - f_7 * pc_z[k] * gsg1_50[k];

        t_171[k] = f_8 * fpf_116[k]
                   + f_4 * pc_x[k] * gpf_116[k];

        t_172[k] = f_1 * gsf_33[k]
                   + f_4 * pc_z[k] * gpf_113[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, pb_z, pc_x, pc_z, fpf_118, fpf_119, \
                         gsg0_55, gsf_36, gsg1_55, gpf_116, gpf_118, \
                         gpf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_8 * fpf_118[k]
                   + f_4 * pc_x[k] * gpf_118[k];

        t_174[k] = f_8 * fpf_119[k]
                   + f_4 * pc_x[k] * gpf_119[k];

        t_175[k] = pb_z[k] * gsg0_55[k]
                   - f_7 * pc_z[k] * gsg1_55[k];

        t_176[k] = f_1 * gsf_36[k]
                   + f_4 * pc_z[k] * gpf_116[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pb_z, pc_y, pc_z, fpf_59, gsg0_57, gsg0_59, \
                         gsf_37, gsf_39, gsg1_57, gsg1_59, gpf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = pb_z[k] * gsg0_57[k]
                   + f_8 * gsf_37[k]
                   - f_7 * pc_z[k] * gsg1_57[k];

        t_178[k] = f_8 * fpf_59[k]
                   + f_4 * pc_y[k] * gpf_119[k];

        t_179[k] = pb_z[k] * gsg0_59[k]
                   + f_0 * gsf_39[k]
                   - f_7 * pc_z[k] * gsg1_59[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pa_y, pa_z, pc_y, pc_z, fpg0_48, fpg0_90, \
                         fpf_30, fpf_60, fpg1_48, fpg1_90, gpf_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = pa_y[k] * fpg0_90[k]
                   - f_7 * pc_y[k] * fpg1_90[k];

        t_181[k] = f_1 * fpf_60[k]
                   + f_4 * pc_y[k] * gpf_120[k];

        t_182[k] = f_1 * fpf_30[k]
                   + f_4 * pc_z[k] * gpf_120[k];

        t_183[k] = pa_z[k] * fpg0_48[k]
                   - f_7 * pc_z[k] * fpg1_48[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, pa_y, pa_z, pc_y, pc_z, fpg0_51, fpg0_95, \
                         fpf_62, fpg1_51, fpg1_95, gpf_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_1 * fpf_62[k]
                   + f_4 * pc_y[k] * gpf_122[k];

        t_185[k] = pa_y[k] * fpg0_95[k]
                   - f_7 * pc_y[k] * fpg1_95[k];

        t_186[k] = pa_z[k] * fpg0_51[k]
                   - f_7 * pc_z[k] * fpg1_51[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, pa_y, pc_x, pc_y, fpg0_99, fpf_127, fpf_128, \
                         fpg1_99, gsf_47, gsf_48, gpf_127, gpf_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = f_8 * fpf_127[k]
                   + f_1 * gsf_47[k]
                   + f_4 * pc_x[k] * gpf_127[k];

        t_188[k] = f_8 * fpf_128[k]
                   + f_1 * gsf_48[k]
                   + f_4 * pc_x[k] * gpf_128[k];

        t_189[k] = pa_y[k] * fpg0_99[k]
                   - f_7 * pc_y[k] * fpg1_99[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, pa_z, pc_y, pc_z, fpg0_55, fpf_36, fpf_68, \
                         fpg1_55, gpd0_77, gpd1_77, gpf_126, gpf_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = pa_z[k] * fpg0_55[k]
                   - f_7 * pc_z[k] * fpg1_55[k];

        t_191[k] = f_1 * fpf_36[k]
                   + f_4 * pc_z[k] * gpf_126[k];

        t_192[k] = f_1 * fpf_68[k]
                   + f_5 * gpd0_77[k]
                   - f_6 * gpd1_77[k]
                   + f_4 * pc_y[k] * gpf_128[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, pa_y, pc_x, pc_y, fpg0_104, fpf_69, fpf_130, \
                         fpg1_104, gpd0_78, gpd1_78, gpf_129, gpf_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_1 * fpf_69[k]
                   + f_4 * pc_y[k] * gpf_129[k];

        t_194[k] = pa_y[k] * fpg0_104[k]
                   - f_7 * pc_y[k] * fpg1_104[k];

        t_195[k] = f_8 * fpf_130[k]
                   + f_2 * gpd0_78[k]
                   - f_3 * gpd1_78[k]
                   + f_4 * pc_x[k] * gpf_130[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pa_z, pc_y, pc_z, fpg0_61, fpg0_63, \
                         fpf_40, fpf_72, fpg1_61, fpg1_63, gsf_42, gpf_130, \
                         gpf_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = pa_z[k] * fpg0_61[k]
                   - f_7 * pc_z[k] * fpg1_61[k];

        t_197[k] = f_1 * fpf_40[k]
                   + f_4 * pc_z[k] * gpf_130[k];

        t_198[k] = pa_z[k] * fpg0_63[k]
                   - f_7 * pc_z[k] * fpg1_63[k];

        t_199[k] = f_1 * fpf_72[k]
                   + f_1 * gsf_42[k]
                   + f_4 * pc_y[k] * gpf_132[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pc_x, fpf_135, fpf_136, fpf_137, fpf_138, \
                         gpd0_83, gpd1_83, gpf_135, gpf_136, gpf_137, \
                         gpf_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_8 * fpf_135[k]
                   + f_5 * gpd0_83[k]
                   - f_6 * gpd1_83[k]
                   + f_4 * pc_x[k] * gpf_135[k];

        t_201[k] = f_8 * fpf_136[k]
                   + f_4 * pc_x[k] * gpf_136[k];

        t_202[k] = f_8 * fpf_137[k]
                   + f_4 * pc_x[k] * gpf_137[k];

        t_203[k] = f_8 * fpf_138[k]
                   + f_4 * pc_x[k] * gpf_138[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, pa_z, pc_x, pc_z, fpg0_70, fpf_46, fpf_139, \
                         fpg1_70, gpf_136, gpf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_8 * fpf_139[k]
                   + f_4 * pc_x[k] * gpf_139[k];

        t_205[k] = pa_z[k] * fpg0_70[k]
                   - f_7 * pc_z[k] * fpg1_70[k];

        t_206[k] = f_1 * fpf_46[k]
                   + f_4 * pc_z[k] * gpf_136[k];
    }

#pragma omp simd aligned(t_207, t_208, pa_x, pc_x, pc_y, dpg0_207, dpg1_207, fpg0_207, fpf_79, \
                         fpg1_207, gsf_49, gpf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_14 * dpg0_207[k]
                   - f_15 * dpg1_207[k]
                   + pa_x[k] * fpg0_207[k]
                   - f_7 * pc_x[k] * fpg1_207[k];

        t_208[k] = f_1 * fpf_79[k]
                   + f_1 * gsf_49[k]
                   + f_4 * pc_y[k] * gpf_139[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, pa_y, pc_y, pc_z, fpg0_120, fpf_49, fpf_80, \
                         fpg1_120, gpd0_83, gpd1_83, gpf_139, gpf_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = f_1 * fpf_49[k]
                   + f_2 * gpd0_83[k]
                   - f_3 * gpd1_83[k]
                   + f_4 * pc_z[k] * gpf_139[k];

        t_210[k] = pa_y[k] * fpg0_120[k]
                   - f_7 * pc_y[k] * fpg1_120[k];

        t_211[k] = f_1 * fpf_80[k]
                   + f_4 * pc_y[k] * gpf_140[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, pa_y, pc_x, pc_y, fpg0_122, fpf_82, fpf_143, \
                         fpg1_122, gpd0_87, gpd1_87, gpf_142, gpf_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = pa_y[k] * fpg0_122[k]
                   - f_7 * pc_y[k] * fpg1_122[k];

        t_213[k] = f_8 * fpf_143[k]
                   + f_5 * gpd0_87[k]
                   - f_6 * gpd1_87[k]
                   + f_4 * pc_x[k] * gpf_143[k];

        t_214[k] = f_1 * fpf_82[k]
                   + f_4 * pc_y[k] * gpf_142[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, pa_y, pc_x, pc_y, fpg0_125, fpf_146, \
                         fpf_147, fpf_148, fpg1_125, gpf_146, gpf_147, \
                         gpf_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = pa_y[k] * fpg0_125[k]
                   - f_7 * pc_y[k] * fpg1_125[k];

        t_216[k] = f_8 * fpf_146[k]
                   + f_4 * pc_x[k] * gpf_146[k];

        t_217[k] = f_8 * fpf_147[k]
                   + f_4 * pc_x[k] * gpf_147[k];

        t_218[k] = f_8 * fpf_148[k]
                   + f_4 * pc_x[k] * gpf_148[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, pc_x, pc_y, pc_z, fpf_56, fpf_86, fpf_149, \
                         gsf_46, gpd0_87, gpd1_87, gpf_146, gpf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_8 * fpf_149[k]
                   + f_4 * pc_x[k] * gpf_149[k];

        t_220[k] = f_1 * fpf_86[k]
                   + f_2 * gpd0_87[k]
                   - f_3 * gpd1_87[k]
                   + f_4 * pc_y[k] * gpf_146[k];

        t_221[k] = f_1 * fpf_56[k]
                   + f_1 * gsf_46[k]
                   + f_4 * pc_z[k] * gpf_146[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, pa_y, pc_y, fpg0_134, fpf_88, fpf_89, fpg1_134, \
                         gpd0_89, gpd1_89, gpf_148, gpf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_1 * fpf_88[k]
                   + f_5 * gpd0_89[k]
                   - f_6 * gpd1_89[k]
                   + f_4 * pc_y[k] * gpf_148[k];

        t_223[k] = f_1 * fpf_89[k]
                   + f_4 * pc_y[k] * gpf_149[k];

        t_224[k] = pa_y[k] * fpg0_134[k]
                   - f_7 * pc_y[k] * fpg1_134[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, t_228, pa_z, pc_y, pc_z, dpg0_0, dpg1_0, \
                         fpg0_90, fpf_60, fpg1_90, gpd0_90, gpd1_90, gpf_150, \
                         gpf_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_14 * dpg0_0[k]
                   - f_15 * dpg1_0[k]
                   + pa_z[k] * fpg0_90[k]
                   - f_7 * pc_z[k] * fpg1_90[k];

        t_226[k] = f_4 * pc_y[k] * gpf_150[k];

        t_227[k] = f_8 * fpf_60[k]
                   + f_4 * pc_z[k] * gpf_150[k];

        t_228[k] = f_5 * gpd0_90[k]
                   - f_6 * gpd1_90[k]
                   + f_4 * pc_y[k] * gpf_151[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, pc_x, pc_y, fpf_155, fpf_156, gsf_55, gsf_56, \
                         gpd0_95, gpd1_95, gpf_152, gpf_155, gpf_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_4 * pc_y[k] * gpf_152[k];

        t_230[k] = f_8 * fpf_155[k]
                   + f_1 * gsf_55[k]
                   + f_5 * gpd0_95[k]
                   - f_6 * gpd1_95[k]
                   + f_4 * pc_x[k] * gpf_155[k];

        t_231[k] = f_8 * fpf_156[k]
                   + f_1 * gsf_56[k]
                   + f_4 * pc_x[k] * gpf_156[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, t_235, pc_x, pc_y, fpf_157, fpf_159, gsf_57, \
                         gsf_59, gpd0_93, gpd1_93, gpf_155, gpf_156, gpf_157, \
                         gpf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_8 * fpf_157[k]
                   + f_1 * gsf_57[k]
                   + f_4 * pc_x[k] * gpf_157[k];

        t_233[k] = f_4 * pc_y[k] * gpf_155[k];

        t_234[k] = f_8 * fpf_159[k]
                   + f_1 * gsf_59[k]
                   + f_4 * pc_x[k] * gpf_159[k];

        t_235[k] = f_2 * gpd0_93[k]
                   - f_3 * gpd1_93[k]
                   + f_4 * pc_y[k] * gpf_156[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, pc_y, pc_z, fpf_69, gpd0_94, gpd0_95, \
                         gpd1_94, gpd1_95, gpf_157, gpf_158, gpf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_12 * gpd0_94[k]
                   - f_13 * gpd1_94[k]
                   + f_4 * pc_y[k] * gpf_157[k];

        t_237[k] = f_5 * gpd0_95[k]
                   - f_6 * gpd1_95[k]
                   + f_4 * pc_y[k] * gpf_158[k];

        t_238[k] = f_4 * pc_y[k] * gpf_159[k];

        t_239[k] = f_8 * fpf_69[k]
                   + f_2 * gpd0_95[k]
                   - f_3 * gpd1_95[k]
                   + f_4 * pc_z[k] * gpf_159[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, pb_y, pc_y, pc_z, fpf_70, gsg0_75, \
                         gsg0_78, gsf_50, gsf_51, gsg1_75, gsg1_78, \
                         gpf_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = pb_y[k] * gsg0_75[k]
                   - f_7 * pc_y[k] * gsg1_75[k];

        t_241[k] = f_1 * gsf_50[k]
                   + f_4 * pc_y[k] * gpf_160[k];

        t_242[k] = f_8 * fpf_70[k]
                   + f_4 * pc_z[k] * gpf_160[k];

        t_243[k] = pb_y[k] * gsg0_78[k]
                   + f_8 * gsf_51[k]
                   - f_7 * pc_y[k] * gsg1_78[k];
    }
}

static auto
compute_prim_gpg_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t dpg0, const size_t dpg1,
                                                          const size_t fpg0, const size_t fpf,
                                                          const size_t fpg1, const size_t gsg0,
                                                          const size_t gsf, const size_t gsg1,
                                                          const size_t gpd0, const size_t gpd1,
                                                          const size_t gpf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 1.5 / gamma;
    const auto f_3 = 1.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = gamma / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 1.5 / q;
    const auto f_12 = 1.0 / gamma;
    const auto f_13 = p / (gamma * q);
    const auto f_14 = 0.5 / p;
    const auto f_15 = 0.5 * gamma / (p * q);

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
    auto *t_362 = buffer.data(target + 362);
    auto *t_363 = buffer.data(target + 363);
    auto *t_364 = buffer.data(target + 364);
    auto *t_365 = buffer.data(target + 365);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dpg0_269 = buffer.data(dpg0 + 269);

    const auto *dpg1_269 = buffer.data(dpg1 + 269);

    const auto *fpg0_135 = buffer.data(fpg0 + 135);
    const auto *fpg0_138 = buffer.data(fpg0 + 138);
    const auto *fpg0_141 = buffer.data(fpg0 + 141);
    const auto *fpg0_150 = buffer.data(fpg0 + 150);
    const auto *fpg0_151 = buffer.data(fpg0 + 151);
    const auto *fpg0_153 = buffer.data(fpg0 + 153);
    const auto *fpg0_225 = buffer.data(fpg0 + 225);
    const auto *fpg0_230 = buffer.data(fpg0 + 230);
    const auto *fpg0_269 = buffer.data(fpg0 + 269);
    const auto *fpg0_285 = buffer.data(fpg0 + 285);
    const auto *fpg0_288 = buffer.data(fpg0 + 288);
    const auto *fpg0_295 = buffer.data(fpg0 + 295);
    const auto *fpg0_297 = buffer.data(fpg0 + 297);
    const auto *fpg0_298 = buffer.data(fpg0 + 298);
    const auto *fpg0_299 = buffer.data(fpg0 + 299);
    const auto *fpg0_305 = buffer.data(fpg0 + 305);
    const auto *fpg0_310 = buffer.data(fpg0 + 310);
    const auto *fpg0_312 = buffer.data(fpg0 + 312);
    const auto *fpg0_314 = buffer.data(fpg0 + 314);
    const auto *fpg0_335 = buffer.data(fpg0 + 335);
    const auto *fpg0_340 = buffer.data(fpg0 + 340);
    const auto *fpg0_342 = buffer.data(fpg0 + 342);
    const auto *fpg0_343 = buffer.data(fpg0 + 343);
    const auto *fpg0_344 = buffer.data(fpg0 + 344);
    const auto *fpg0_350 = buffer.data(fpg0 + 350);
    const auto *fpg0_355 = buffer.data(fpg0 + 355);
    const auto *fpg0_356 = buffer.data(fpg0 + 356);
    const auto *fpg0_357 = buffer.data(fpg0 + 357);
    const auto *fpg0_359 = buffer.data(fpg0 + 359);

    const auto *fpf_80 = buffer.data(fpf + 80);
    const auto *fpf_90 = buffer.data(fpf + 90);
    const auto *fpf_96 = buffer.data(fpf + 96);
    const auto *fpf_99 = buffer.data(fpf + 99);
    const auto *fpf_100 = buffer.data(fpf + 100);
    const auto *fpf_106 = buffer.data(fpf + 106);
    const auto *fpf_110 = buffer.data(fpf + 110);
    const auto *fpf_119 = buffer.data(fpf + 119);
    const auto *fpf_120 = buffer.data(fpf + 120);
    const auto *fpf_122 = buffer.data(fpf + 122);
    const auto *fpf_126 = buffer.data(fpf + 126);
    const auto *fpf_128 = buffer.data(fpf + 128);
    const auto *fpf_129 = buffer.data(fpf + 129);
    const auto *fpf_132 = buffer.data(fpf + 132);
    const auto *fpf_140 = buffer.data(fpf + 140);
    const auto *fpf_142 = buffer.data(fpf + 142);
    const auto *fpf_149 = buffer.data(fpf + 149);
    const auto *fpf_150 = buffer.data(fpf + 150);
    const auto *fpf_152 = buffer.data(fpf + 152);
    const auto *fpf_166 = buffer.data(fpf + 166);
    const auto *fpf_167 = buffer.data(fpf + 167);
    const auto *fpf_169 = buffer.data(fpf + 169);
    const auto *fpf_170 = buffer.data(fpf + 170);
    const auto *fpf_175 = buffer.data(fpf + 175);
    const auto *fpf_176 = buffer.data(fpf + 176);
    const auto *fpf_177 = buffer.data(fpf + 177);
    const auto *fpf_179 = buffer.data(fpf + 179);
    const auto *fpf_180 = buffer.data(fpf + 180);
    const auto *fpf_183 = buffer.data(fpf + 183);
    const auto *fpf_186 = buffer.data(fpf + 186);
    const auto *fpf_188 = buffer.data(fpf + 188);
    const auto *fpf_189 = buffer.data(fpf + 189);
    const auto *fpf_190 = buffer.data(fpf + 190);
    const auto *fpf_193 = buffer.data(fpf + 193);
    const auto *fpf_196 = buffer.data(fpf + 196);
    const auto *fpf_198 = buffer.data(fpf + 198);
    const auto *fpf_199 = buffer.data(fpf + 199);
    const auto *fpf_205 = buffer.data(fpf + 205);
    const auto *fpf_206 = buffer.data(fpf + 206);
    const auto *fpf_208 = buffer.data(fpf + 208);
    const auto *fpf_209 = buffer.data(fpf + 209);
    const auto *fpf_215 = buffer.data(fpf + 215);
    const auto *fpf_217 = buffer.data(fpf + 217);
    const auto *fpf_218 = buffer.data(fpf + 218);
    const auto *fpf_219 = buffer.data(fpf + 219);
    const auto *fpf_225 = buffer.data(fpf + 225);
    const auto *fpf_226 = buffer.data(fpf + 226);
    const auto *fpf_227 = buffer.data(fpf + 227);
    const auto *fpf_228 = buffer.data(fpf + 228);
    const auto *fpf_229 = buffer.data(fpf + 229);
    const auto *fpf_230 = buffer.data(fpf + 230);
    const auto *fpf_233 = buffer.data(fpf + 233);
    const auto *fpf_235 = buffer.data(fpf + 235);
    const auto *fpf_236 = buffer.data(fpf + 236);
    const auto *fpf_237 = buffer.data(fpf + 237);
    const auto *fpf_238 = buffer.data(fpf + 238);
    const auto *fpf_239 = buffer.data(fpf + 239);
    const auto *fpf_243 = buffer.data(fpf + 243);

    const auto *fpg1_135 = buffer.data(fpg1 + 135);
    const auto *fpg1_138 = buffer.data(fpg1 + 138);
    const auto *fpg1_141 = buffer.data(fpg1 + 141);
    const auto *fpg1_150 = buffer.data(fpg1 + 150);
    const auto *fpg1_151 = buffer.data(fpg1 + 151);
    const auto *fpg1_153 = buffer.data(fpg1 + 153);
    const auto *fpg1_225 = buffer.data(fpg1 + 225);
    const auto *fpg1_230 = buffer.data(fpg1 + 230);
    const auto *fpg1_269 = buffer.data(fpg1 + 269);
    const auto *fpg1_285 = buffer.data(fpg1 + 285);
    const auto *fpg1_288 = buffer.data(fpg1 + 288);
    const auto *fpg1_295 = buffer.data(fpg1 + 295);
    const auto *fpg1_297 = buffer.data(fpg1 + 297);
    const auto *fpg1_298 = buffer.data(fpg1 + 298);
    const auto *fpg1_299 = buffer.data(fpg1 + 299);
    const auto *fpg1_305 = buffer.data(fpg1 + 305);
    const auto *fpg1_310 = buffer.data(fpg1 + 310);
    const auto *fpg1_312 = buffer.data(fpg1 + 312);
    const auto *fpg1_314 = buffer.data(fpg1 + 314);
    const auto *fpg1_335 = buffer.data(fpg1 + 335);
    const auto *fpg1_340 = buffer.data(fpg1 + 340);
    const auto *fpg1_342 = buffer.data(fpg1 + 342);
    const auto *fpg1_343 = buffer.data(fpg1 + 343);
    const auto *fpg1_344 = buffer.data(fpg1 + 344);
    const auto *fpg1_350 = buffer.data(fpg1 + 350);
    const auto *fpg1_355 = buffer.data(fpg1 + 355);
    const auto *fpg1_356 = buffer.data(fpg1 + 356);
    const auto *fpg1_357 = buffer.data(fpg1 + 357);
    const auto *fpg1_359 = buffer.data(fpg1 + 359);

    const auto *gsg0_80 = buffer.data(gsg0 + 80);
    const auto *gsg0_85 = buffer.data(gsg0 + 85);
    const auto *gsg0_86 = buffer.data(gsg0 + 86);
    const auto *gsg0_87 = buffer.data(gsg0 + 87);
    const auto *gsg0_89 = buffer.data(gsg0 + 89);
    const auto *gsg0_90 = buffer.data(gsg0 + 90);
    const auto *gsg0_93 = buffer.data(gsg0 + 93);

    const auto *gsf_50 = buffer.data(gsf + 50);
    const auto *gsf_52 = buffer.data(gsf + 52);
    const auto *gsf_55 = buffer.data(gsf + 55);
    const auto *gsf_56 = buffer.data(gsf + 56);
    const auto *gsf_57 = buffer.data(gsf + 57);
    const auto *gsf_58 = buffer.data(gsf + 58);
    const auto *gsf_59 = buffer.data(gsf + 59);
    const auto *gsf_60 = buffer.data(gsf + 60);
    const auto *gsf_61 = buffer.data(gsf + 61);
    const auto *gsf_63 = buffer.data(gsf + 63);
    const auto *gsf_66 = buffer.data(gsf + 66);
    const auto *gsf_68 = buffer.data(gsf + 68);
    const auto *gsf_69 = buffer.data(gsf + 69);
    const auto *gsf_70 = buffer.data(gsf + 70);
    const auto *gsf_72 = buffer.data(gsf + 72);
    const auto *gsf_75 = buffer.data(gsf + 75);
    const auto *gsf_77 = buffer.data(gsf + 77);
    const auto *gsf_78 = buffer.data(gsf + 78);
    const auto *gsf_79 = buffer.data(gsf + 79);
    const auto *gsf_83 = buffer.data(gsf + 83);

    const auto *gsg1_80 = buffer.data(gsg1 + 80);
    const auto *gsg1_85 = buffer.data(gsg1 + 85);
    const auto *gsg1_86 = buffer.data(gsg1 + 86);
    const auto *gsg1_87 = buffer.data(gsg1 + 87);
    const auto *gsg1_89 = buffer.data(gsg1 + 89);
    const auto *gsg1_90 = buffer.data(gsg1 + 90);
    const auto *gsg1_93 = buffer.data(gsg1 + 93);

    const auto *gpd0_102 = buffer.data(gpd0 + 102);
    const auto *gpd0_105 = buffer.data(gpd0 + 105);
    const auto *gpd0_106 = buffer.data(gpd0 + 106);
    const auto *gpd0_107 = buffer.data(gpd0 + 107);
    const auto *gpd0_108 = buffer.data(gpd0 + 108);
    const auto *gpd0_111 = buffer.data(gpd0 + 111);
    const auto *gpd0_113 = buffer.data(gpd0 + 113);
    const auto *gpd0_114 = buffer.data(gpd0 + 114);
    const auto *gpd0_129 = buffer.data(gpd0 + 129);
    const auto *gpd0_131 = buffer.data(gpd0 + 131);
    const auto *gpd0_138 = buffer.data(gpd0 + 138);
    const auto *gpd0_141 = buffer.data(gpd0 + 141);
    const auto *gpd0_147 = buffer.data(gpd0 + 147);

    const auto *gpd1_102 = buffer.data(gpd1 + 102);
    const auto *gpd1_105 = buffer.data(gpd1 + 105);
    const auto *gpd1_106 = buffer.data(gpd1 + 106);
    const auto *gpd1_107 = buffer.data(gpd1 + 107);
    const auto *gpd1_108 = buffer.data(gpd1 + 108);
    const auto *gpd1_111 = buffer.data(gpd1 + 111);
    const auto *gpd1_113 = buffer.data(gpd1 + 113);
    const auto *gpd1_114 = buffer.data(gpd1 + 114);
    const auto *gpd1_129 = buffer.data(gpd1 + 129);
    const auto *gpd1_131 = buffer.data(gpd1 + 131);
    const auto *gpd1_138 = buffer.data(gpd1 + 138);
    const auto *gpd1_141 = buffer.data(gpd1 + 141);
    const auto *gpd1_147 = buffer.data(gpd1 + 147);

    const auto *gpf_162 = buffer.data(gpf + 162);
    const auto *gpf_165 = buffer.data(gpf + 165);
    const auto *gpf_166 = buffer.data(gpf + 166);
    const auto *gpf_167 = buffer.data(gpf + 167);
    const auto *gpf_169 = buffer.data(gpf + 169);
    const auto *gpf_170 = buffer.data(gpf + 170);
    const auto *gpf_171 = buffer.data(gpf + 171);
    const auto *gpf_172 = buffer.data(gpf + 172);
    const auto *gpf_175 = buffer.data(gpf + 175);
    const auto *gpf_176 = buffer.data(gpf + 176);
    const auto *gpf_177 = buffer.data(gpf + 177);
    const auto *gpf_178 = buffer.data(gpf + 178);
    const auto *gpf_179 = buffer.data(gpf + 179);
    const auto *gpf_180 = buffer.data(gpf + 180);
    const auto *gpf_181 = buffer.data(gpf + 181);
    const auto *gpf_182 = buffer.data(gpf + 182);
    const auto *gpf_183 = buffer.data(gpf + 183);
    const auto *gpf_186 = buffer.data(gpf + 186);
    const auto *gpf_187 = buffer.data(gpf + 187);
    const auto *gpf_188 = buffer.data(gpf + 188);
    const auto *gpf_189 = buffer.data(gpf + 189);
    const auto *gpf_190 = buffer.data(gpf + 190);
    const auto *gpf_191 = buffer.data(gpf + 191);
    const auto *gpf_192 = buffer.data(gpf + 192);
    const auto *gpf_193 = buffer.data(gpf + 193);
    const auto *gpf_196 = buffer.data(gpf + 196);
    const auto *gpf_198 = buffer.data(gpf + 198);
    const auto *gpf_199 = buffer.data(gpf + 199);
    const auto *gpf_200 = buffer.data(gpf + 200);
    const auto *gpf_201 = buffer.data(gpf + 201);
    const auto *gpf_203 = buffer.data(gpf + 203);
    const auto *gpf_206 = buffer.data(gpf + 206);
    const auto *gpf_208 = buffer.data(gpf + 208);
    const auto *gpf_209 = buffer.data(gpf + 209);
    const auto *gpf_210 = buffer.data(gpf + 210);
    const auto *gpf_212 = buffer.data(gpf + 212);
    const auto *gpf_215 = buffer.data(gpf + 215);
    const auto *gpf_216 = buffer.data(gpf + 216);
    const auto *gpf_217 = buffer.data(gpf + 217);
    const auto *gpf_218 = buffer.data(gpf + 218);
    const auto *gpf_219 = buffer.data(gpf + 219);
    const auto *gpf_220 = buffer.data(gpf + 220);
    const auto *gpf_222 = buffer.data(gpf + 222);
    const auto *gpf_226 = buffer.data(gpf + 226);
    const auto *gpf_227 = buffer.data(gpf + 227);
    const auto *gpf_228 = buffer.data(gpf + 228);
    const auto *gpf_229 = buffer.data(gpf + 229);
    const auto *gpf_230 = buffer.data(gpf + 230);
    const auto *gpf_232 = buffer.data(gpf + 232);
    const auto *gpf_233 = buffer.data(gpf + 233);
    const auto *gpf_236 = buffer.data(gpf + 236);
    const auto *gpf_237 = buffer.data(gpf + 237);
    const auto *gpf_238 = buffer.data(gpf + 238);
    const auto *gpf_239 = buffer.data(gpf + 239);
    const auto *gpf_240 = buffer.data(gpf + 240);
    const auto *gpf_242 = buffer.data(gpf + 242);
    const auto *gpf_243 = buffer.data(gpf + 243);

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pb_y, pc_x, pc_y, fpf_166, fpf_167, \
                         gsg0_80, gsf_52, gsg1_80, gpf_162, gpf_166, \
                         gpf_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_1 * gsf_52[k]
                   + f_4 * pc_y[k] * gpf_162[k];

        t_245[k] = pb_y[k] * gsg0_80[k]
                   - f_7 * pc_y[k] * gsg1_80[k];

        t_246[k] = f_8 * fpf_166[k]
                   + f_4 * pc_x[k] * gpf_166[k];

        t_247[k] = f_8 * fpf_167[k]
                   + f_4 * pc_x[k] * gpf_167[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, pb_y, pc_x, pc_y, fpf_169, gsg0_85, gsf_55, \
                         gsf_56, gsg1_85, gpf_165, gpf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_1 * gsf_55[k]
                   + f_4 * pc_y[k] * gpf_165[k];

        t_249[k] = f_8 * fpf_169[k]
                   + f_4 * pc_x[k] * gpf_169[k];

        t_250[k] = pb_y[k] * gsg0_85[k]
                   + f_0 * gsf_56[k]
                   - f_7 * pc_y[k] * gsg1_85[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pb_y, pc_y, gsg0_86, gsg0_87, gsg0_89, \
                         gsf_57, gsf_58, gsf_59, gsg1_86, gsg1_87, gsg1_89, \
                         gpf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = pb_y[k] * gsg0_86[k]
                   + f_9 * gsf_57[k]
                   - f_7 * pc_y[k] * gsg1_86[k];

        t_252[k] = pb_y[k] * gsg0_87[k]
                   + f_8 * gsf_58[k]
                   - f_7 * pc_y[k] * gsg1_87[k];

        t_253[k] = f_1 * gsf_59[k]
                   + f_4 * pc_y[k] * gpf_169[k];

        t_254[k] = pb_y[k] * gsg0_89[k]
                   - f_7 * pc_y[k] * gsg1_89[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, pc_x, pc_y, pc_z, fpf_80, fpf_170, \
                         gsf_50, gpd0_102, gpd1_102, gpf_170, gpf_171, \
                         gpf_172 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_8 * fpf_170[k]
                   + f_2 * gpd0_102[k]
                   - f_3 * gpd1_102[k]
                   + f_4 * pc_x[k] * gpf_170[k];

        t_256[k] = f_4 * pc_y[k] * gpf_170[k];

        t_257[k] = f_8 * fpf_80[k]
                   + f_1 * gsf_50[k]
                   + f_4 * pc_z[k] * gpf_170[k];

        t_258[k] = f_5 * gpd0_102[k]
                   - f_6 * gpd1_102[k]
                   + f_4 * pc_y[k] * gpf_171[k];

        t_259[k] = f_4 * pc_y[k] * gpf_172[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pc_x, pc_y, fpf_175, fpf_176, fpf_177, \
                         gpd0_107, gpd1_107, gpf_175, gpf_176, \
                         gpf_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_8 * fpf_175[k]
                   + f_5 * gpd0_107[k]
                   - f_6 * gpd1_107[k]
                   + f_4 * pc_x[k] * gpf_175[k];

        t_261[k] = f_8 * fpf_176[k]
                   + f_4 * pc_x[k] * gpf_176[k];

        t_262[k] = f_8 * fpf_177[k]
                   + f_4 * pc_x[k] * gpf_177[k];

        t_263[k] = f_4 * pc_y[k] * gpf_175[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pc_x, pc_y, fpf_179, gpd0_105, gpd0_106, \
                         gpd1_105, gpd1_106, gpf_176, gpf_177, \
                         gpf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_8 * fpf_179[k]
                   + f_4 * pc_x[k] * gpf_179[k];

        t_265[k] = f_2 * gpd0_105[k]
                   - f_3 * gpd1_105[k]
                   + f_4 * pc_y[k] * gpf_176[k];

        t_266[k] = f_12 * gpd0_106[k]
                   - f_13 * gpd1_106[k]
                   + f_4 * pc_y[k] * gpf_177[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, pa_x, pc_x, pc_y, dpg0_269, dpg1_269, fpg0_269, \
                         fpg1_269, gpd0_107, gpd1_107, gpf_178, \
                         gpf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_5 * gpd0_107[k]
                   - f_6 * gpd1_107[k]
                   + f_4 * pc_y[k] * gpf_178[k];

        t_268[k] = f_4 * pc_y[k] * gpf_179[k];

        t_269[k] = f_14 * dpg0_269[k]
                   - f_15 * dpg1_269[k]
                   + pa_x[k] * fpg0_269[k]
                   - f_7 * pc_x[k] * fpg1_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, pc_x, pc_y, pc_z, fpf_90, fpf_180, gsf_60, \
                         gpd0_108, gpd1_108, gpf_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_1 * fpf_180[k]
                   + f_1 * gsf_60[k]
                   + f_2 * gpd0_108[k]
                   - f_3 * gpd1_108[k]
                   + f_4 * pc_x[k] * gpf_180[k];

        t_271[k] = f_9 * fpf_90[k]
                   + f_4 * pc_y[k] * gpf_180[k];

        t_272[k] = f_4 * pc_z[k] * gpf_180[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, pc_x, pc_z, fpf_183, gsf_63, gpd0_108, gpd0_111, \
                         gpd1_108, gpd1_111, gpf_181, gpf_182, \
                         gpf_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_1 * fpf_183[k]
                   + f_1 * gsf_63[k]
                   + f_5 * gpd0_111[k]
                   - f_6 * gpd1_111[k]
                   + f_4 * pc_x[k] * gpf_183[k];

        t_274[k] = f_4 * pc_z[k] * gpf_181[k];

        t_275[k] = f_5 * gpd0_108[k]
                   - f_6 * gpd1_108[k]
                   + f_4 * pc_z[k] * gpf_182[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_x, pc_z, fpf_186, fpf_188, fpf_189, \
                         gsf_66, gsf_68, gsf_69, gpf_183, gpf_186, gpf_188, \
                         gpf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_1 * fpf_186[k]
                   + f_1 * gsf_66[k]
                   + f_4 * pc_x[k] * gpf_186[k];

        t_277[k] = f_4 * pc_z[k] * gpf_183[k];

        t_278[k] = f_1 * fpf_188[k]
                   + f_1 * gsf_68[k]
                   + f_4 * pc_x[k] * gpf_188[k];

        t_279[k] = f_1 * fpf_189[k]
                   + f_1 * gsf_69[k]
                   + f_4 * pc_x[k] * gpf_189[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, pc_y, pc_z, fpf_96, fpf_99, \
                         gpd0_111, gpd0_113, gpd1_111, gpd1_113, gpf_186, gpf_187, \
                         gpf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_9 * fpf_96[k]
                   + f_2 * gpd0_111[k]
                   - f_3 * gpd1_111[k]
                   + f_4 * pc_y[k] * gpf_186[k];

        t_281[k] = f_4 * pc_z[k] * gpf_186[k];

        t_282[k] = f_5 * gpd0_111[k]
                   - f_6 * gpd1_111[k]
                   + f_4 * pc_z[k] * gpf_187[k];

        t_283[k] = f_9 * fpf_99[k]
                   + f_4 * pc_y[k] * gpf_189[k];

        t_284[k] = f_2 * gpd0_113[k]
                   - f_3 * gpd1_113[k]
                   + f_4 * pc_z[k] * gpf_189[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, pa_x, pc_x, pc_y, pc_z, fpg0_285, fpf_100, \
                         fpf_190, fpg1_285, gsf_60, gpf_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = pa_x[k] * fpg0_285[k]
                   + f_0 * fpf_190[k]
                   - f_7 * pc_x[k] * fpg1_285[k];

        t_286[k] = f_9 * fpf_100[k]
                   + f_1 * gsf_60[k]
                   + f_4 * pc_y[k] * gpf_190[k];

        t_287[k] = f_4 * pc_z[k] * gpf_190[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, t_291, pa_x, pc_x, pc_z, fpg0_288, fpf_193, \
                         fpf_196, fpg1_288, gpd0_114, gpd1_114, gpf_191, gpf_192, \
                         gpf_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = pa_x[k] * fpg0_288[k]
                   + f_8 * fpf_193[k]
                   - f_7 * pc_x[k] * fpg1_288[k];

        t_289[k] = f_4 * pc_z[k] * gpf_191[k];

        t_290[k] = f_5 * gpd0_114[k]
                   - f_6 * gpd1_114[k]
                   + f_4 * pc_z[k] * gpf_192[k];

        t_291[k] = f_1 * fpf_196[k]
                   + f_4 * pc_x[k] * gpf_196[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, t_296, pa_x, pc_x, pc_z, fpg0_295, \
                         fpf_198, fpf_199, fpg1_295, gpf_193, gpf_196, gpf_198, \
                         gpf_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = f_4 * pc_z[k] * gpf_193[k];

        t_293[k] = f_1 * fpf_198[k]
                   + f_4 * pc_x[k] * gpf_198[k];

        t_294[k] = f_1 * fpf_199[k]
                   + f_4 * pc_x[k] * gpf_199[k];

        t_295[k] = pa_x[k] * fpg0_295[k]
                   - f_7 * pc_x[k] * fpg1_295[k];

        t_296[k] = f_4 * pc_z[k] * gpf_196[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, t_300, pa_x, pb_z, pc_x, pc_z, fpg0_297, \
                         fpg0_298, fpg0_299, fpg1_297, fpg1_298, fpg1_299, gsg0_90, \
                         gsg1_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = pa_x[k] * fpg0_297[k]
                   - f_7 * pc_x[k] * fpg1_297[k];

        t_298[k] = pa_x[k] * fpg0_298[k]
                   - f_7 * pc_x[k] * fpg1_298[k];

        t_299[k] = pa_x[k] * fpg0_299[k]
                   - f_7 * pc_x[k] * fpg1_299[k];

        t_300[k] = pb_z[k] * gsg0_90[k]
                   - f_7 * pc_z[k] * gsg1_90[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pb_z, pc_y, pc_z, fpf_110, gsg0_93, \
                         gsf_60, gsf_61, gsg1_93, gpf_200, gpf_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_9 * fpf_110[k]
                   + f_4 * pc_y[k] * gpf_200[k];

        t_302[k] = f_1 * gsf_60[k]
                   + f_4 * pc_z[k] * gpf_200[k];

        t_303[k] = pb_z[k] * gsg0_93[k]
                   - f_7 * pc_z[k] * gsg1_93[k];

        t_304[k] = f_1 * gsf_61[k]
                   + f_4 * pc_z[k] * gpf_201[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pa_x, pc_x, pc_z, fpg0_305, fpf_205, \
                         fpf_206, fpf_208, fpg1_305, gsf_63, gpf_203, gpf_206, \
                         gpf_208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = pa_x[k] * fpg0_305[k]
                   + f_8 * fpf_205[k]
                   - f_7 * pc_x[k] * fpg1_305[k];

        t_306[k] = f_1 * fpf_206[k]
                   + f_4 * pc_x[k] * gpf_206[k];

        t_307[k] = f_1 * gsf_63[k]
                   + f_4 * pc_z[k] * gpf_203[k];

        t_308[k] = f_1 * fpf_208[k]
                   + f_4 * pc_x[k] * gpf_208[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, pa_x, pc_x, pc_z, fpg0_310, fpg0_312, \
                         fpf_209, fpg1_310, fpg1_312, gsf_66, gpf_206, \
                         gpf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_1 * fpf_209[k]
                   + f_4 * pc_x[k] * gpf_209[k];

        t_310[k] = pa_x[k] * fpg0_310[k]
                   - f_7 * pc_x[k] * fpg1_310[k];

        t_311[k] = f_1 * gsf_66[k]
                   + f_4 * pc_z[k] * gpf_206[k];

        t_312[k] = pa_x[k] * fpg0_312[k]
                   - f_7 * pc_x[k] * fpg1_312[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, pa_x, pa_z, pc_x, pc_y, pc_z, fpg0_135, \
                         fpg0_314, fpf_119, fpg1_135, fpg1_314, \
                         gpf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = f_9 * fpf_119[k]
                   + f_4 * pc_y[k] * gpf_209[k];

        t_314[k] = pa_x[k] * fpg0_314[k]
                   - f_7 * pc_x[k] * fpg1_314[k];

        t_315[k] = pa_z[k] * fpg0_135[k]
                   - f_7 * pc_z[k] * fpg1_135[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, t_319, pa_z, pc_y, pc_z, fpg0_138, fpf_90, \
                         fpf_120, fpf_122, fpg1_138, gpf_210, gpf_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = f_8 * fpf_120[k]
                   + f_4 * pc_y[k] * gpf_210[k];

        t_317[k] = f_1 * fpf_90[k]
                   + f_4 * pc_z[k] * gpf_210[k];

        t_318[k] = pa_z[k] * fpg0_138[k]
                   - f_7 * pc_z[k] * fpg1_138[k];

        t_319[k] = f_8 * fpf_122[k]
                   + f_4 * pc_y[k] * gpf_212[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, pa_z, pc_x, pc_z, fpg0_141, fpf_215, fpf_217, \
                         fpg1_141, gsf_75, gsf_77, gpd0_131, gpd1_131, gpf_215, \
                         gpf_217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_1 * fpf_215[k]
                   + f_1 * gsf_75[k]
                   + f_5 * gpd0_131[k]
                   - f_6 * gpd1_131[k]
                   + f_4 * pc_x[k] * gpf_215[k];

        t_321[k] = pa_z[k] * fpg0_141[k]
                   - f_7 * pc_z[k] * fpg1_141[k];

        t_322[k] = f_1 * fpf_217[k]
                   + f_1 * gsf_77[k]
                   + f_4 * pc_x[k] * gpf_217[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, pc_x, pc_y, fpf_126, fpf_218, fpf_219, gsf_78, \
                         gsf_79, gpd0_129, gpd1_129, gpf_216, gpf_218, \
                         gpf_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = f_1 * fpf_218[k]
                   + f_1 * gsf_78[k]
                   + f_4 * pc_x[k] * gpf_218[k];

        t_324[k] = f_1 * fpf_219[k]
                   + f_1 * gsf_79[k]
                   + f_4 * pc_x[k] * gpf_219[k];

        t_325[k] = f_8 * fpf_126[k]
                   + f_2 * gpd0_129[k]
                   - f_3 * gpd1_129[k]
                   + f_4 * pc_y[k] * gpf_216[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, pc_y, pc_z, fpf_96, fpf_99, fpf_128, \
                         fpf_129, gpd0_131, gpd1_131, gpf_216, gpf_218, \
                         gpf_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_1 * fpf_96[k]
                   + f_4 * pc_z[k] * gpf_216[k];

        t_327[k] = f_8 * fpf_128[k]
                   + f_5 * gpd0_131[k]
                   - f_6 * gpd1_131[k]
                   + f_4 * pc_y[k] * gpf_218[k];

        t_328[k] = f_8 * fpf_129[k]
                   + f_4 * pc_y[k] * gpf_219[k];

        t_329[k] = f_1 * fpf_99[k]
                   + f_2 * gpd0_131[k]
                   - f_3 * gpd1_131[k]
                   + f_4 * pc_z[k] * gpf_219[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, pa_z, pc_z, fpg0_150, fpg0_151, fpg0_153, \
                         fpf_100, fpg1_150, fpg1_151, fpg1_153, \
                         gpf_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = pa_z[k] * fpg0_150[k]
                   - f_7 * pc_z[k] * fpg1_150[k];

        t_331[k] = pa_z[k] * fpg0_151[k]
                   - f_7 * pc_z[k] * fpg1_151[k];

        t_332[k] = f_1 * fpf_100[k]
                   + f_4 * pc_z[k] * gpf_220[k];

        t_333[k] = pa_z[k] * fpg0_153[k]
                   - f_7 * pc_z[k] * fpg1_153[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, pa_x, pc_x, pc_y, fpg0_335, fpf_132, fpf_225, \
                         fpf_226, fpg1_335, gsf_72, gpf_222, gpf_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_8 * fpf_132[k]
                   + f_1 * gsf_72[k]
                   + f_4 * pc_y[k] * gpf_222[k];

        t_335[k] = pa_x[k] * fpg0_335[k]
                   + f_8 * fpf_225[k]
                   - f_7 * pc_x[k] * fpg1_335[k];

        t_336[k] = f_1 * fpf_226[k]
                   + f_4 * pc_x[k] * gpf_226[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, t_340, pa_x, pc_x, fpg0_340, fpf_227, fpf_228, \
                         fpf_229, fpg1_340, gpf_227, gpf_228, gpf_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = f_1 * fpf_227[k]
                   + f_4 * pc_x[k] * gpf_227[k];

        t_338[k] = f_1 * fpf_228[k]
                   + f_4 * pc_x[k] * gpf_228[k];

        t_339[k] = f_1 * fpf_229[k]
                   + f_4 * pc_x[k] * gpf_229[k];

        t_340[k] = pa_x[k] * fpg0_340[k]
                   - f_7 * pc_x[k] * fpg1_340[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, t_344, pa_x, pc_x, pc_z, fpg0_342, fpg0_343, \
                         fpg0_344, fpf_106, fpg1_342, fpg1_343, fpg1_344, \
                         gpf_226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_1 * fpf_106[k]
                   + f_4 * pc_z[k] * gpf_226[k];

        t_342[k] = pa_x[k] * fpg0_342[k]
                   - f_7 * pc_x[k] * fpg1_342[k];

        t_343[k] = pa_x[k] * fpg0_343[k]
                   - f_7 * pc_x[k] * fpg1_343[k];

        t_344[k] = pa_x[k] * fpg0_344[k]
                   - f_7 * pc_x[k] * fpg1_344[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, pc_x, pc_y, pc_z, fpf_110, fpf_140, fpf_230, \
                         gsf_70, gpd0_138, gpd1_138, gpf_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = f_1 * fpf_230[k]
                   + f_2 * gpd0_138[k]
                   - f_3 * gpd1_138[k]
                   + f_4 * pc_x[k] * gpf_230[k];

        t_346[k] = f_8 * fpf_140[k]
                   + f_4 * pc_y[k] * gpf_230[k];

        t_347[k] = f_1 * fpf_110[k]
                   + f_1 * gsf_70[k]
                   + f_4 * pc_z[k] * gpf_230[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, pa_x, pc_x, pc_y, fpg0_350, fpf_142, fpf_233, \
                         fpf_235, fpg1_350, gpd0_141, gpd1_141, gpf_232, \
                         gpf_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_1 * fpf_233[k]
                   + f_5 * gpd0_141[k]
                   - f_6 * gpd1_141[k]
                   + f_4 * pc_x[k] * gpf_233[k];

        t_349[k] = f_8 * fpf_142[k]
                   + f_4 * pc_y[k] * gpf_232[k];

        t_350[k] = pa_x[k] * fpg0_350[k]
                   + f_8 * fpf_235[k]
                   - f_7 * pc_x[k] * fpg1_350[k];
    }

#pragma omp simd aligned(t_351, t_352, t_353, t_354, pc_x, fpf_236, fpf_237, fpf_238, fpf_239, \
                         gpf_236, gpf_237, gpf_238, gpf_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_351[k] = f_1 * fpf_236[k]
                   + f_4 * pc_x[k] * gpf_236[k];

        t_352[k] = f_1 * fpf_237[k]
                   + f_4 * pc_x[k] * gpf_237[k];

        t_353[k] = f_1 * fpf_238[k]
                   + f_4 * pc_x[k] * gpf_238[k];

        t_354[k] = f_1 * fpf_239[k]
                   + f_4 * pc_x[k] * gpf_239[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, pa_x, pc_x, pc_y, fpg0_355, fpg0_356, \
                         fpg0_357, fpf_149, fpg1_355, fpg1_356, fpg1_357, \
                         gpf_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = pa_x[k] * fpg0_355[k]
                   - f_7 * pc_x[k] * fpg1_355[k];

        t_356[k] = pa_x[k] * fpg0_356[k]
                   - f_7 * pc_x[k] * fpg1_356[k];

        t_357[k] = pa_x[k] * fpg0_357[k]
                   - f_7 * pc_x[k] * fpg1_357[k];

        t_358[k] = f_8 * fpf_149[k]
                   + f_4 * pc_y[k] * gpf_239[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, t_362, pa_x, pa_y, pc_x, pc_y, pc_z, fpg0_225, \
                         fpg0_359, fpf_120, fpf_150, fpg1_225, fpg1_359, \
                         gpf_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = pa_x[k] * fpg0_359[k]
                   - f_7 * pc_x[k] * fpg1_359[k];

        t_360[k] = pa_y[k] * fpg0_225[k]
                   - f_7 * pc_y[k] * fpg1_225[k];

        t_361[k] = f_1 * fpf_150[k]
                   + f_4 * pc_y[k] * gpf_240[k];

        t_362[k] = f_8 * fpf_120[k]
                   + f_4 * pc_z[k] * gpf_240[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, pa_y, pc_x, pc_y, fpg0_230, fpf_152, fpf_243, \
                         fpg1_230, gsf_83, gpd0_147, gpd1_147, gpf_242, \
                         gpf_243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = f_1 * fpf_243[k]
                   + f_1 * gsf_83[k]
                   + f_5 * gpd0_147[k]
                   - f_6 * gpd1_147[k]
                   + f_4 * pc_x[k] * gpf_243[k];

        t_364[k] = f_1 * fpf_152[k]
                   + f_4 * pc_y[k] * gpf_242[k];

        t_365[k] = pa_y[k] * fpg0_230[k]
                   - f_7 * pc_y[k] * fpg1_230[k];
    }
}

static auto
compute_prim_gpg_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t fpg0, const size_t fpf,
                                                          const size_t fpg1, const size_t gsg0,
                                                          const size_t gsf, const size_t gsg1,
                                                          const size_t gpd0, const size_t gpd1,
                                                          const size_t gpf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 1.5 / gamma;
    const auto f_3 = 1.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = gamma / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 1.5 / q;
    const auto f_12 = 1.0 / gamma;
    const auto f_13 = p / (gamma * q);

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
    auto *t_488 = buffer.data(target + 488);
    auto *t_489 = buffer.data(target + 489);
    auto *t_490 = buffer.data(target + 490);
    auto *t_491 = buffer.data(target + 491);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *fpg0_234 = buffer.data(fpg0 + 234);
    const auto *fpg0_255 = buffer.data(fpg0 + 255);
    const auto *fpg0_257 = buffer.data(fpg0 + 257);
    const auto *fpg0_260 = buffer.data(fpg0 + 260);
    const auto *fpg0_378 = buffer.data(fpg0 + 378);
    const auto *fpg0_385 = buffer.data(fpg0 + 385);
    const auto *fpg0_387 = buffer.data(fpg0 + 387);
    const auto *fpg0_388 = buffer.data(fpg0 + 388);
    const auto *fpg0_389 = buffer.data(fpg0 + 389);
    const auto *fpg0_393 = buffer.data(fpg0 + 393);
    const auto *fpg0_400 = buffer.data(fpg0 + 400);
    const auto *fpg0_401 = buffer.data(fpg0 + 401);
    const auto *fpg0_402 = buffer.data(fpg0 + 402);
    const auto *fpg0_404 = buffer.data(fpg0 + 404);
    const auto *fpg0_423 = buffer.data(fpg0 + 423);
    const auto *fpg0_430 = buffer.data(fpg0 + 430);
    const auto *fpg0_431 = buffer.data(fpg0 + 431);
    const auto *fpg0_432 = buffer.data(fpg0 + 432);
    const auto *fpg0_434 = buffer.data(fpg0 + 434);
    const auto *fpg0_435 = buffer.data(fpg0 + 435);
    const auto *fpg0_440 = buffer.data(fpg0 + 440);
    const auto *fpg0_445 = buffer.data(fpg0 + 445);
    const auto *fpg0_446 = buffer.data(fpg0 + 446);
    const auto *fpg0_447 = buffer.data(fpg0 + 447);
    const auto *fpg0_449 = buffer.data(fpg0 + 449);

    const auto *fpf_126 = buffer.data(fpf + 126);
    const auto *fpf_129 = buffer.data(fpf + 129);
    const auto *fpf_130 = buffer.data(fpf + 130);
    const auto *fpf_136 = buffer.data(fpf + 136);
    const auto *fpf_150 = buffer.data(fpf + 150);
    const auto *fpf_156 = buffer.data(fpf + 156);
    const auto *fpf_158 = buffer.data(fpf + 158);
    const auto *fpf_159 = buffer.data(fpf + 159);
    const auto *fpf_160 = buffer.data(fpf + 160);
    const auto *fpf_162 = buffer.data(fpf + 162);
    const auto *fpf_170 = buffer.data(fpf + 170);
    const auto *fpf_172 = buffer.data(fpf + 172);
    const auto *fpf_179 = buffer.data(fpf + 179);
    const auto *fpf_189 = buffer.data(fpf + 189);
    const auto *fpf_196 = buffer.data(fpf + 196);
    const auto *fpf_199 = buffer.data(fpf + 199);
    const auto *fpf_246 = buffer.data(fpf + 246);
    const auto *fpf_247 = buffer.data(fpf + 247);
    const auto *fpf_248 = buffer.data(fpf + 248);
    const auto *fpf_250 = buffer.data(fpf + 250);
    const auto *fpf_253 = buffer.data(fpf + 253);
    const auto *fpf_255 = buffer.data(fpf + 255);
    const auto *fpf_256 = buffer.data(fpf + 256);
    const auto *fpf_257 = buffer.data(fpf + 257);
    const auto *fpf_258 = buffer.data(fpf + 258);
    const auto *fpf_259 = buffer.data(fpf + 259);
    const auto *fpf_263 = buffer.data(fpf + 263);
    const auto *fpf_266 = buffer.data(fpf + 266);
    const auto *fpf_267 = buffer.data(fpf + 267);
    const auto *fpf_268 = buffer.data(fpf + 268);
    const auto *fpf_269 = buffer.data(fpf + 269);
    const auto *fpf_270 = buffer.data(fpf + 270);
    const auto *fpf_275 = buffer.data(fpf + 275);
    const auto *fpf_276 = buffer.data(fpf + 276);
    const auto *fpf_277 = buffer.data(fpf + 277);
    const auto *fpf_279 = buffer.data(fpf + 279);
    const auto *fpf_283 = buffer.data(fpf + 283);
    const auto *fpf_286 = buffer.data(fpf + 286);
    const auto *fpf_287 = buffer.data(fpf + 287);
    const auto *fpf_289 = buffer.data(fpf + 289);
    const auto *fpf_290 = buffer.data(fpf + 290);
    const auto *fpf_295 = buffer.data(fpf + 295);
    const auto *fpf_296 = buffer.data(fpf + 296);
    const auto *fpf_297 = buffer.data(fpf + 297);
    const auto *fpf_299 = buffer.data(fpf + 299);

    const auto *fpg1_234 = buffer.data(fpg1 + 234);
    const auto *fpg1_255 = buffer.data(fpg1 + 255);
    const auto *fpg1_257 = buffer.data(fpg1 + 257);
    const auto *fpg1_260 = buffer.data(fpg1 + 260);
    const auto *fpg1_378 = buffer.data(fpg1 + 378);
    const auto *fpg1_385 = buffer.data(fpg1 + 385);
    const auto *fpg1_387 = buffer.data(fpg1 + 387);
    const auto *fpg1_388 = buffer.data(fpg1 + 388);
    const auto *fpg1_389 = buffer.data(fpg1 + 389);
    const auto *fpg1_393 = buffer.data(fpg1 + 393);
    const auto *fpg1_400 = buffer.data(fpg1 + 400);
    const auto *fpg1_401 = buffer.data(fpg1 + 401);
    const auto *fpg1_402 = buffer.data(fpg1 + 402);
    const auto *fpg1_404 = buffer.data(fpg1 + 404);
    const auto *fpg1_423 = buffer.data(fpg1 + 423);
    const auto *fpg1_430 = buffer.data(fpg1 + 430);
    const auto *fpg1_431 = buffer.data(fpg1 + 431);
    const auto *fpg1_432 = buffer.data(fpg1 + 432);
    const auto *fpg1_434 = buffer.data(fpg1 + 434);
    const auto *fpg1_435 = buffer.data(fpg1 + 435);
    const auto *fpg1_440 = buffer.data(fpg1 + 440);
    const auto *fpg1_445 = buffer.data(fpg1 + 445);
    const auto *fpg1_446 = buffer.data(fpg1 + 446);
    const auto *fpg1_447 = buffer.data(fpg1 + 447);
    const auto *fpg1_449 = buffer.data(fpg1 + 449);

    const auto *gsg0_135 = buffer.data(gsg0 + 135);
    const auto *gsg0_140 = buffer.data(gsg0 + 140);
    const auto *gsg0_150 = buffer.data(gsg0 + 150);
    const auto *gsg0_151 = buffer.data(gsg0 + 151);
    const auto *gsg0_153 = buffer.data(gsg0 + 153);
    const auto *gsg0_155 = buffer.data(gsg0 + 155);
    const auto *gsg0_160 = buffer.data(gsg0 + 160);
    const auto *gsg0_162 = buffer.data(gsg0 + 162);
    const auto *gsg0_164 = buffer.data(gsg0 + 164);

    const auto *gsf_80 = buffer.data(gsf + 80);
    const auto *gsf_82 = buffer.data(gsf + 82);
    const auto *gsf_86 = buffer.data(gsf + 86);
    const auto *gsf_87 = buffer.data(gsf + 87);
    const auto *gsf_88 = buffer.data(gsf + 88);
    const auto *gsf_90 = buffer.data(gsf + 90);
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

    const auto *gsg1_135 = buffer.data(gsg1 + 135);
    const auto *gsg1_140 = buffer.data(gsg1 + 140);
    const auto *gsg1_150 = buffer.data(gsg1 + 150);
    const auto *gsg1_151 = buffer.data(gsg1 + 151);
    const auto *gsg1_153 = buffer.data(gsg1 + 153);
    const auto *gsg1_155 = buffer.data(gsg1 + 155);
    const auto *gsg1_160 = buffer.data(gsg1 + 160);
    const auto *gsg1_162 = buffer.data(gsg1 + 162);
    const auto *gsg1_164 = buffer.data(gsg1 + 164);

    const auto *gpd0_147 = buffer.data(gpd0 + 147);
    const auto *gpd0_149 = buffer.data(gpd0 + 149);
    const auto *gpd0_150 = buffer.data(gpd0 + 150);
    const auto *gpd0_155 = buffer.data(gpd0 + 155);
    const auto *gpd0_162 = buffer.data(gpd0 + 162);
    const auto *gpd0_165 = buffer.data(gpd0 + 165);
    const auto *gpd0_166 = buffer.data(gpd0 + 166);
    const auto *gpd0_167 = buffer.data(gpd0 + 167);
    const auto *gpd0_174 = buffer.data(gpd0 + 174);
    const auto *gpd0_186 = buffer.data(gpd0 + 186);
    const auto *gpd0_187 = buffer.data(gpd0 + 187);
    const auto *gpd0_189 = buffer.data(gpd0 + 189);
    const auto *gpd0_191 = buffer.data(gpd0 + 191);
    const auto *gpd0_197 = buffer.data(gpd0 + 197);

    const auto *gpd1_147 = buffer.data(gpd1 + 147);
    const auto *gpd1_149 = buffer.data(gpd1 + 149);
    const auto *gpd1_150 = buffer.data(gpd1 + 150);
    const auto *gpd1_155 = buffer.data(gpd1 + 155);
    const auto *gpd1_162 = buffer.data(gpd1 + 162);
    const auto *gpd1_165 = buffer.data(gpd1 + 165);
    const auto *gpd1_166 = buffer.data(gpd1 + 166);
    const auto *gpd1_167 = buffer.data(gpd1 + 167);
    const auto *gpd1_174 = buffer.data(gpd1 + 174);
    const auto *gpd1_186 = buffer.data(gpd1 + 186);
    const auto *gpd1_187 = buffer.data(gpd1 + 187);
    const auto *gpd1_189 = buffer.data(gpd1 + 189);
    const auto *gpd1_191 = buffer.data(gpd1 + 191);
    const auto *gpd1_197 = buffer.data(gpd1 + 197);

    const auto *gpf_246 = buffer.data(gpf + 246);
    const auto *gpf_247 = buffer.data(gpf + 247);
    const auto *gpf_248 = buffer.data(gpf + 248);
    const auto *gpf_249 = buffer.data(gpf + 249);
    const auto *gpf_250 = buffer.data(gpf + 250);
    const auto *gpf_252 = buffer.data(gpf + 252);
    const auto *gpf_255 = buffer.data(gpf + 255);
    const auto *gpf_256 = buffer.data(gpf + 256);
    const auto *gpf_257 = buffer.data(gpf + 257);
    const auto *gpf_258 = buffer.data(gpf + 258);
    const auto *gpf_259 = buffer.data(gpf + 259);
    const auto *gpf_260 = buffer.data(gpf + 260);
    const auto *gpf_262 = buffer.data(gpf + 262);
    const auto *gpf_266 = buffer.data(gpf + 266);
    const auto *gpf_267 = buffer.data(gpf + 267);
    const auto *gpf_268 = buffer.data(gpf + 268);
    const auto *gpf_269 = buffer.data(gpf + 269);
    const auto *gpf_270 = buffer.data(gpf + 270);
    const auto *gpf_271 = buffer.data(gpf + 271);
    const auto *gpf_272 = buffer.data(gpf + 272);
    const auto *gpf_275 = buffer.data(gpf + 275);
    const auto *gpf_276 = buffer.data(gpf + 276);
    const auto *gpf_277 = buffer.data(gpf + 277);
    const auto *gpf_278 = buffer.data(gpf + 278);
    const auto *gpf_279 = buffer.data(gpf + 279);
    const auto *gpf_280 = buffer.data(gpf + 280);
    const auto *gpf_282 = buffer.data(gpf + 282);
    const auto *gpf_285 = buffer.data(gpf + 285);
    const auto *gpf_286 = buffer.data(gpf + 286);
    const auto *gpf_287 = buffer.data(gpf + 287);
    const auto *gpf_289 = buffer.data(gpf + 289);
    const auto *gpf_290 = buffer.data(gpf + 290);
    const auto *gpf_291 = buffer.data(gpf + 291);
    const auto *gpf_292 = buffer.data(gpf + 292);
    const auto *gpf_295 = buffer.data(gpf + 295);
    const auto *gpf_296 = buffer.data(gpf + 296);
    const auto *gpf_297 = buffer.data(gpf + 297);
    const auto *gpf_299 = buffer.data(gpf + 299);
    const auto *gpf_300 = buffer.data(gpf + 300);
    const auto *gpf_301 = buffer.data(gpf + 301);
    const auto *gpf_306 = buffer.data(gpf + 306);
    const auto *gpf_307 = buffer.data(gpf + 307);
    const auto *gpf_308 = buffer.data(gpf + 308);
    const auto *gpf_309 = buffer.data(gpf + 309);
    const auto *gpf_310 = buffer.data(gpf + 310);
    const auto *gpf_311 = buffer.data(gpf + 311);
    const auto *gpf_313 = buffer.data(gpf + 313);
    const auto *gpf_315 = buffer.data(gpf + 315);
    const auto *gpf_316 = buffer.data(gpf + 316);
    const auto *gpf_317 = buffer.data(gpf + 317);
    const auto *gpf_318 = buffer.data(gpf + 318);
    const auto *gpf_319 = buffer.data(gpf + 319);
    const auto *gpf_320 = buffer.data(gpf + 320);
    const auto *gpf_321 = buffer.data(gpf + 321);
    const auto *gpf_325 = buffer.data(gpf + 325);
    const auto *gpf_326 = buffer.data(gpf + 326);
    const auto *gpf_327 = buffer.data(gpf + 327);
    const auto *gpf_328 = buffer.data(gpf + 328);
    const auto *gpf_329 = buffer.data(gpf + 329);

#pragma omp simd aligned(t_366, t_367, t_368, pc_x, fpf_246, fpf_247, fpf_248, gsf_86, gsf_87, \
                         gsf_88, gpf_246, gpf_247, gpf_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_1 * fpf_246[k]
                   + f_1 * gsf_86[k]
                   + f_4 * pc_x[k] * gpf_246[k];

        t_367[k] = f_1 * fpf_247[k]
                   + f_1 * gsf_87[k]
                   + f_4 * pc_x[k] * gpf_247[k];

        t_368[k] = f_1 * fpf_248[k]
                   + f_1 * gsf_88[k]
                   + f_4 * pc_x[k] * gpf_248[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, pa_y, pc_y, pc_z, fpg0_234, fpf_126, fpf_156, \
                         fpg1_234, gpd0_147, gpd1_147, gpf_246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = pa_y[k] * fpg0_234[k]
                   - f_7 * pc_y[k] * fpg1_234[k];

        t_370[k] = f_1 * fpf_156[k]
                   + f_2 * gpd0_147[k]
                   - f_3 * gpd1_147[k]
                   + f_4 * pc_y[k] * gpf_246[k];

        t_371[k] = f_8 * fpf_126[k]
                   + f_4 * pc_z[k] * gpf_246[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, pc_y, pc_z, fpf_129, fpf_158, fpf_159, gpd0_149, \
                         gpd1_149, gpf_248, gpf_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = f_1 * fpf_158[k]
                   + f_5 * gpd0_149[k]
                   - f_6 * gpd1_149[k]
                   + f_4 * pc_y[k] * gpf_248[k];

        t_373[k] = f_1 * fpf_159[k]
                   + f_4 * pc_y[k] * gpf_249[k];

        t_374[k] = f_8 * fpf_129[k]
                   + f_2 * gpd0_149[k]
                   - f_3 * gpd1_149[k]
                   + f_4 * pc_z[k] * gpf_249[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, pc_x, pc_y, pc_z, fpf_130, fpf_160, fpf_250, \
                         gsf_80, gpd0_150, gpd1_150, gpf_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = f_1 * fpf_250[k]
                   + f_2 * gpd0_150[k]
                   - f_3 * gpd1_150[k]
                   + f_4 * pc_x[k] * gpf_250[k];

        t_376[k] = f_1 * fpf_160[k]
                   + f_1 * gsf_80[k]
                   + f_4 * pc_y[k] * gpf_250[k];

        t_377[k] = f_8 * fpf_130[k]
                   + f_4 * pc_z[k] * gpf_250[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, pa_x, pc_x, pc_y, fpg0_378, fpf_162, fpf_253, \
                         fpf_255, fpg1_378, gsf_82, gpd0_155, gpd1_155, gpf_252, \
                         gpf_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = pa_x[k] * fpg0_378[k]
                   + f_8 * fpf_253[k]
                   - f_7 * pc_x[k] * fpg1_378[k];

        t_379[k] = f_1 * fpf_162[k]
                   + f_1 * gsf_82[k]
                   + f_4 * pc_y[k] * gpf_252[k];

        t_380[k] = f_1 * fpf_255[k]
                   + f_5 * gpd0_155[k]
                   - f_6 * gpd1_155[k]
                   + f_4 * pc_x[k] * gpf_255[k];
    }

#pragma omp simd aligned(t_381, t_382, t_383, t_384, pc_x, fpf_256, fpf_257, fpf_258, fpf_259, \
                         gpf_256, gpf_257, gpf_258, gpf_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_381[k] = f_1 * fpf_256[k]
                   + f_4 * pc_x[k] * gpf_256[k];

        t_382[k] = f_1 * fpf_257[k]
                   + f_4 * pc_x[k] * gpf_257[k];

        t_383[k] = f_1 * fpf_258[k]
                   + f_4 * pc_x[k] * gpf_258[k];

        t_384[k] = f_1 * fpf_259[k]
                   + f_4 * pc_x[k] * gpf_259[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, t_388, pa_x, pc_x, pc_z, fpg0_385, fpg0_387, \
                         fpg0_388, fpf_136, fpg1_385, fpg1_387, fpg1_388, \
                         gpf_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = pa_x[k] * fpg0_385[k]
                   - f_7 * pc_x[k] * fpg1_385[k];

        t_386[k] = f_8 * fpf_136[k]
                   + f_4 * pc_z[k] * gpf_256[k];

        t_387[k] = pa_x[k] * fpg0_387[k]
                   - f_7 * pc_x[k] * fpg1_387[k];

        t_388[k] = pa_x[k] * fpg0_388[k]
                   - f_7 * pc_x[k] * fpg1_388[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, t_392, pa_x, pa_y, pc_x, pc_y, fpg0_255, \
                         fpg0_257, fpg0_389, fpf_170, fpg1_255, fpg1_257, fpg1_389, \
                         gpf_260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = pa_x[k] * fpg0_389[k]
                   - f_7 * pc_x[k] * fpg1_389[k];

        t_390[k] = pa_y[k] * fpg0_255[k]
                   - f_7 * pc_y[k] * fpg1_255[k];

        t_391[k] = f_1 * fpf_170[k]
                   + f_4 * pc_y[k] * gpf_260[k];

        t_392[k] = pa_y[k] * fpg0_257[k]
                   - f_7 * pc_y[k] * fpg1_257[k];
    }

#pragma omp simd aligned(t_393, t_394, t_395, pa_x, pa_y, pc_x, pc_y, fpg0_260, fpg0_393, \
                         fpf_172, fpf_263, fpg1_260, fpg1_393, \
                         gpf_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_393[k] = pa_x[k] * fpg0_393[k]
                   + f_8 * fpf_263[k]
                   - f_7 * pc_x[k] * fpg1_393[k];

        t_394[k] = f_1 * fpf_172[k]
                   + f_4 * pc_y[k] * gpf_262[k];

        t_395[k] = pa_y[k] * fpg0_260[k]
                   - f_7 * pc_y[k] * fpg1_260[k];
    }

#pragma omp simd aligned(t_396, t_397, t_398, t_399, pc_x, fpf_266, fpf_267, fpf_268, fpf_269, \
                         gpf_266, gpf_267, gpf_268, gpf_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_396[k] = f_1 * fpf_266[k]
                   + f_4 * pc_x[k] * gpf_266[k];

        t_397[k] = f_1 * fpf_267[k]
                   + f_4 * pc_x[k] * gpf_267[k];

        t_398[k] = f_1 * fpf_268[k]
                   + f_4 * pc_x[k] * gpf_268[k];

        t_399[k] = f_1 * fpf_269[k]
                   + f_4 * pc_x[k] * gpf_269[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, pa_x, pc_x, pc_y, fpg0_400, fpg0_401, \
                         fpg0_402, fpf_179, fpg1_400, fpg1_401, fpg1_402, \
                         gpf_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = pa_x[k] * fpg0_400[k]
                   - f_7 * pc_x[k] * fpg1_400[k];

        t_401[k] = pa_x[k] * fpg0_401[k]
                   - f_7 * pc_x[k] * fpg1_401[k];

        t_402[k] = pa_x[k] * fpg0_402[k]
                   - f_7 * pc_x[k] * fpg1_402[k];

        t_403[k] = f_1 * fpf_179[k]
                   + f_4 * pc_y[k] * gpf_269[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, pa_x, pc_x, pc_y, pc_z, fpg0_404, \
                         fpf_150, fpf_270, fpg1_404, gsf_90, gpd0_162, gpd1_162, \
                         gpf_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = pa_x[k] * fpg0_404[k]
                   - f_7 * pc_x[k] * fpg1_404[k];

        t_405[k] = f_1 * fpf_270[k]
                   + f_1 * gsf_90[k]
                   + f_2 * gpd0_162[k]
                   - f_3 * gpd1_162[k]
                   + f_4 * pc_x[k] * gpf_270[k];

        t_406[k] = f_4 * pc_y[k] * gpf_270[k];

        t_407[k] = f_9 * fpf_150[k]
                   + f_4 * pc_z[k] * gpf_270[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, pc_x, pc_y, fpf_275, gsf_95, gpd0_162, gpd0_167, \
                         gpd1_162, gpd1_167, gpf_271, gpf_272, \
                         gpf_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = f_5 * gpd0_162[k]
                   - f_6 * gpd1_162[k]
                   + f_4 * pc_y[k] * gpf_271[k];

        t_409[k] = f_4 * pc_y[k] * gpf_272[k];

        t_410[k] = f_1 * fpf_275[k]
                   + f_1 * gsf_95[k]
                   + f_5 * gpd0_167[k]
                   - f_6 * gpd1_167[k]
                   + f_4 * pc_x[k] * gpf_275[k];
    }

#pragma omp simd aligned(t_411, t_412, t_413, t_414, pc_x, pc_y, fpf_276, fpf_277, fpf_279, \
                         gsf_96, gsf_97, gsf_99, gpf_275, gpf_276, gpf_277, \
                         gpf_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_411[k] = f_1 * fpf_276[k]
                   + f_1 * gsf_96[k]
                   + f_4 * pc_x[k] * gpf_276[k];

        t_412[k] = f_1 * fpf_277[k]
                   + f_1 * gsf_97[k]
                   + f_4 * pc_x[k] * gpf_277[k];

        t_413[k] = f_4 * pc_y[k] * gpf_275[k];

        t_414[k] = f_1 * fpf_279[k]
                   + f_1 * gsf_99[k]
                   + f_4 * pc_x[k] * gpf_279[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, t_418, pc_y, gpd0_165, gpd0_166, gpd0_167, \
                         gpd1_165, gpd1_166, gpd1_167, gpf_276, gpf_277, gpf_278, \
                         gpf_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = f_2 * gpd0_165[k]
                   - f_3 * gpd1_165[k]
                   + f_4 * pc_y[k] * gpf_276[k];

        t_416[k] = f_12 * gpd0_166[k]
                   - f_13 * gpd1_166[k]
                   + f_4 * pc_y[k] * gpf_277[k];

        t_417[k] = f_5 * gpd0_167[k]
                   - f_6 * gpd1_167[k]
                   + f_4 * pc_y[k] * gpf_278[k];

        t_418[k] = f_4 * pc_y[k] * gpf_279[k];
    }

#pragma omp simd aligned(t_419, t_420, t_421, t_422, pb_y, pc_y, pc_z, fpf_159, fpf_160, \
                         gsg0_135, gsf_90, gsg1_135, gpd0_167, gpd1_167, gpf_279, \
                         gpf_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_419[k] = f_9 * fpf_159[k]
                   + f_2 * gpd0_167[k]
                   - f_3 * gpd1_167[k]
                   + f_4 * pc_z[k] * gpf_279[k];

        t_420[k] = pb_y[k] * gsg0_135[k]
                   - f_7 * pc_y[k] * gsg1_135[k];

        t_421[k] = f_1 * gsf_90[k]
                   + f_4 * pc_y[k] * gpf_280[k];

        t_422[k] = f_9 * fpf_160[k]
                   + f_4 * pc_z[k] * gpf_280[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, pa_x, pb_y, pc_x, pc_y, fpg0_423, fpf_283, \
                         fpg1_423, gsg0_140, gsf_92, gsg1_140, \
                         gpf_282 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = pa_x[k] * fpg0_423[k]
                   + f_8 * fpf_283[k]
                   - f_7 * pc_x[k] * fpg1_423[k];

        t_424[k] = f_1 * gsf_92[k]
                   + f_4 * pc_y[k] * gpf_282[k];

        t_425[k] = pb_y[k] * gsg0_140[k]
                   - f_7 * pc_y[k] * gsg1_140[k];
    }

#pragma omp simd aligned(t_426, t_427, t_428, t_429, pc_x, pc_y, fpf_286, fpf_287, fpf_289, \
                         gsf_95, gpf_285, gpf_286, gpf_287, gpf_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_426[k] = f_1 * fpf_286[k]
                   + f_4 * pc_x[k] * gpf_286[k];

        t_427[k] = f_1 * fpf_287[k]
                   + f_4 * pc_x[k] * gpf_287[k];

        t_428[k] = f_1 * gsf_95[k]
                   + f_4 * pc_y[k] * gpf_285[k];

        t_429[k] = f_1 * fpf_289[k]
                   + f_4 * pc_x[k] * gpf_289[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, pa_x, pc_x, pc_y, fpg0_430, fpg0_431, \
                         fpg0_432, fpg1_430, fpg1_431, fpg1_432, gsf_99, \
                         gpf_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = pa_x[k] * fpg0_430[k]
                   - f_7 * pc_x[k] * fpg1_430[k];

        t_431[k] = pa_x[k] * fpg0_431[k]
                   - f_7 * pc_x[k] * fpg1_431[k];

        t_432[k] = pa_x[k] * fpg0_432[k]
                   - f_7 * pc_x[k] * fpg1_432[k];

        t_433[k] = f_1 * gsf_99[k]
                   + f_4 * pc_y[k] * gpf_289[k];
    }

#pragma omp simd aligned(t_434, t_435, t_436, t_437, pa_x, pc_x, pc_y, pc_z, fpg0_434, \
                         fpg0_435, fpf_170, fpf_290, fpg1_434, fpg1_435, gsf_90, \
                         gpf_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_434[k] = pa_x[k] * fpg0_434[k]
                   - f_7 * pc_x[k] * fpg1_434[k];

        t_435[k] = pa_x[k] * fpg0_435[k]
                   + f_0 * fpf_290[k]
                   - f_7 * pc_x[k] * fpg1_435[k];

        t_436[k] = f_4 * pc_y[k] * gpf_290[k];

        t_437[k] = f_9 * fpf_170[k]
                   + f_1 * gsf_90[k]
                   + f_4 * pc_z[k] * gpf_290[k];
    }

#pragma omp simd aligned(t_438, t_439, t_440, t_441, pa_x, pc_x, pc_y, fpg0_440, fpf_295, \
                         fpf_296, fpg1_440, gpd0_174, gpd1_174, gpf_291, gpf_292, \
                         gpf_296 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_438[k] = f_5 * gpd0_174[k]
                   - f_6 * gpd1_174[k]
                   + f_4 * pc_y[k] * gpf_291[k];

        t_439[k] = f_4 * pc_y[k] * gpf_292[k];

        t_440[k] = pa_x[k] * fpg0_440[k]
                   + f_8 * fpf_295[k]
                   - f_7 * pc_x[k] * fpg1_440[k];

        t_441[k] = f_1 * fpf_296[k]
                   + f_4 * pc_x[k] * gpf_296[k];
    }

#pragma omp simd aligned(t_442, t_443, t_444, t_445, pa_x, pc_x, pc_y, fpg0_445, fpf_297, \
                         fpf_299, fpg1_445, gpf_295, gpf_297, gpf_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_442[k] = f_1 * fpf_297[k]
                   + f_4 * pc_x[k] * gpf_297[k];

        t_443[k] = f_4 * pc_y[k] * gpf_295[k];

        t_444[k] = f_1 * fpf_299[k]
                   + f_4 * pc_x[k] * gpf_299[k];

        t_445[k] = pa_x[k] * fpg0_445[k]
                   - f_7 * pc_x[k] * fpg1_445[k];
    }

#pragma omp simd aligned(t_446, t_447, t_448, t_449, pa_x, pc_x, pc_y, fpg0_446, fpg0_447, \
                         fpg0_449, fpg1_446, fpg1_447, fpg1_449, \
                         gpf_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_446[k] = pa_x[k] * fpg0_446[k]
                   - f_7 * pc_x[k] * fpg1_446[k];

        t_447[k] = pa_x[k] * fpg0_447[k]
                   - f_7 * pc_x[k] * fpg1_447[k];

        t_448[k] = f_4 * pc_y[k] * gpf_299[k];

        t_449[k] = pa_x[k] * fpg0_449[k]
                   - f_7 * pc_x[k] * fpg1_449[k];
    }

#pragma omp simd aligned(t_450, t_451, t_452, pb_x, pc_x, pc_z, gsg0_150, gsg0_151, gsf_100, \
                         gsf_101, gsg1_150, gsg1_151, gpf_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_450[k] = pb_x[k] * gsg0_150[k]
                   + f_0 * gsf_100[k]
                   - f_7 * pc_x[k] * gsg1_150[k];

        t_451[k] = pb_x[k] * gsg0_151[k]
                   + f_9 * gsf_101[k]
                   - f_7 * pc_x[k] * gsg1_151[k];

        t_452[k] = f_4 * pc_z[k] * gpf_300[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, t_456, pb_x, pc_x, pc_z, gsg0_153, gsg0_155, \
                         gsf_103, gsf_105, gsf_106, gsg1_153, gsg1_155, gpf_301, \
                         gpf_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = pb_x[k] * gsg0_153[k]
                   + f_8 * gsf_103[k]
                   - f_7 * pc_x[k] * gsg1_153[k];

        t_454[k] = f_4 * pc_z[k] * gpf_301[k];

        t_455[k] = pb_x[k] * gsg0_155[k]
                   + f_8 * gsf_105[k]
                   - f_7 * pc_x[k] * gsg1_155[k];

        t_456[k] = f_1 * gsf_106[k]
                   + f_4 * pc_x[k] * gpf_306[k];
    }

#pragma omp simd aligned(t_457, t_458, t_459, t_460, pb_x, pc_x, gsg0_160, gsf_107, gsf_108, \
                         gsf_109, gsg1_160, gpf_307, gpf_308, gpf_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_457[k] = f_1 * gsf_107[k]
                   + f_4 * pc_x[k] * gpf_307[k];

        t_458[k] = f_1 * gsf_108[k]
                   + f_4 * pc_x[k] * gpf_308[k];

        t_459[k] = f_1 * gsf_109[k]
                   + f_4 * pc_x[k] * gpf_309[k];

        t_460[k] = pb_x[k] * gsg0_160[k]
                   - f_7 * pc_x[k] * gsg1_160[k];
    }

#pragma omp simd aligned(t_461, t_462, t_463, t_464, pb_x, pc_x, pc_y, pc_z, fpf_189, \
                         gsg0_162, gsg0_164, gsg1_162, gsg1_164, gpf_306, \
                         gpf_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_461[k] = f_4 * pc_z[k] * gpf_306[k];

        t_462[k] = pb_x[k] * gsg0_162[k]
                   - f_7 * pc_x[k] * gsg1_162[k];

        t_463[k] = f_0 * fpf_189[k]
                   + f_4 * pc_y[k] * gpf_309[k];

        t_464[k] = pb_x[k] * gsg0_164[k]
                   - f_7 * pc_x[k] * gsg1_164[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, t_469, pc_x, pc_z, gpd0_186, gpd0_187, \
                         gpd0_189, gpd1_186, gpd1_187, gpd1_189, gpf_310, gpf_311, \
                         gpf_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = f_2 * gpd0_186[k]
                   - f_3 * gpd1_186[k]
                   + f_4 * pc_x[k] * gpf_310[k];

        t_466[k] = f_12 * gpd0_187[k]
                   - f_13 * gpd1_187[k]
                   + f_4 * pc_x[k] * gpf_311[k];

        t_467[k] = f_4 * pc_z[k] * gpf_310[k];

        t_468[k] = f_5 * gpd0_189[k]
                   - f_6 * gpd1_189[k]
                   + f_4 * pc_x[k] * gpf_313[k];

        t_469[k] = f_4 * pc_z[k] * gpf_311[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, t_474, pc_x, gpd0_191, gpd1_191, gpf_315, \
                         gpf_316, gpf_317, gpf_318, gpf_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_5 * gpd0_191[k]
                   - f_6 * gpd1_191[k]
                   + f_4 * pc_x[k] * gpf_315[k];

        t_471[k] = f_4 * pc_x[k] * gpf_316[k];

        t_472[k] = f_4 * pc_x[k] * gpf_317[k];

        t_473[k] = f_4 * pc_x[k] * gpf_318[k];

        t_474[k] = f_4 * pc_x[k] * gpf_319[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, pc_y, pc_z, fpf_196, fpf_199, gsf_106, \
                         gsf_109, gpd0_189, gpd1_189, gpf_316, gpf_317, \
                         gpf_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = f_0 * fpf_196[k]
                   + f_1 * gsf_106[k]
                   + f_2 * gpd0_189[k]
                   - f_3 * gpd1_189[k]
                   + f_4 * pc_y[k] * gpf_316[k];

        t_476[k] = f_4 * pc_z[k] * gpf_316[k];

        t_477[k] = f_5 * gpd0_189[k]
                   - f_6 * gpd1_189[k]
                   + f_4 * pc_z[k] * gpf_317[k];

        t_478[k] = f_0 * fpf_199[k]
                   + f_1 * gsf_109[k]
                   + f_4 * pc_y[k] * gpf_319[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, t_482, pb_z, pc_z, gsg0_150, gsg0_151, gsf_100, \
                         gsg1_150, gsg1_151, gpd0_191, gpd1_191, gpf_319, \
                         gpf_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = f_2 * gpd0_191[k]
                   - f_3 * gpd1_191[k]
                   + f_4 * pc_z[k] * gpf_319[k];

        t_480[k] = pb_z[k] * gsg0_150[k]
                   - f_7 * pc_z[k] * gsg1_150[k];

        t_481[k] = pb_z[k] * gsg0_151[k]
                   - f_7 * pc_z[k] * gsg1_151[k];

        t_482[k] = f_1 * gsf_100[k]
                   + f_4 * pc_z[k] * gpf_320[k];
    }

#pragma omp simd aligned(t_483, t_484, t_485, t_486, pb_z, pc_x, pc_z, gsg0_153, gsf_101, \
                         gsg1_153, gpd0_197, gpd1_197, gpf_321, gpf_325, \
                         gpf_326 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_483[k] = pb_z[k] * gsg0_153[k]
                   - f_7 * pc_z[k] * gsg1_153[k];

        t_484[k] = f_1 * gsf_101[k]
                   + f_4 * pc_z[k] * gpf_321[k];

        t_485[k] = f_5 * gpd0_197[k]
                   - f_6 * gpd1_197[k]
                   + f_4 * pc_x[k] * gpf_325[k];

        t_486[k] = f_4 * pc_x[k] * gpf_326[k];
    }

#pragma omp simd aligned(t_487, t_488, t_489, t_490, t_491, pb_z, pc_x, pc_z, gsg0_160, \
                         gsf_106, gsg1_160, gpf_326, gpf_327, gpf_328, \
                         gpf_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_487[k] = f_4 * pc_x[k] * gpf_327[k];

        t_488[k] = f_4 * pc_x[k] * gpf_328[k];

        t_489[k] = f_4 * pc_x[k] * gpf_329[k];

        t_490[k] = pb_z[k] * gsg0_160[k]
                   - f_7 * pc_z[k] * gsg1_160[k];

        t_491[k] = f_1 * gsf_106[k]
                   + f_4 * pc_z[k] * gpf_326[k];
    }
}

static auto
compute_prim_gpg_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t dpg0, const size_t dpg1,
                                                          const size_t fpg0, const size_t fpf,
                                                          const size_t fpg1, const size_t gsg0,
                                                          const size_t gsf, const size_t gsg1,
                                                          const size_t gpd0, const size_t gpd1,
                                                          const size_t gpf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 1.5 / gamma;
    const auto f_3 = 1.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = gamma / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 1.5 / q;
    const auto f_10 = 1.0 / p;
    const auto f_11 = gamma / (p * q);
    const auto f_12 = 1.0 / gamma;
    const auto f_13 = p / (gamma * q);
    const auto f_14 = 0.5 / p;
    const auto f_15 = 0.5 * gamma / (p * q);

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dpg0_160 = buffer.data(dpg0 + 160);
    const auto *dpg0_224 = buffer.data(dpg0 + 224);
    const auto *dpg0_269 = buffer.data(dpg0 + 269);

    const auto *dpg1_160 = buffer.data(dpg1 + 160);
    const auto *dpg1_224 = buffer.data(dpg1 + 224);
    const auto *dpg1_269 = buffer.data(dpg1 + 269);

    const auto *fpg0_270 = buffer.data(fpg0 + 270);
    const auto *fpg0_271 = buffer.data(fpg0 + 271);
    const auto *fpg0_273 = buffer.data(fpg0 + 273);
    const auto *fpg0_274 = buffer.data(fpg0 + 274);
    const auto *fpg0_280 = buffer.data(fpg0 + 280);
    const auto *fpg0_285 = buffer.data(fpg0 + 285);
    const auto *fpg0_286 = buffer.data(fpg0 + 286);
    const auto *fpg0_288 = buffer.data(fpg0 + 288);
    const auto *fpg0_295 = buffer.data(fpg0 + 295);
    const auto *fpg0_297 = buffer.data(fpg0 + 297);
    const auto *fpg0_340 = buffer.data(fpg0 + 340);
    const auto *fpg0_359 = buffer.data(fpg0 + 359);
    const auto *fpg0_404 = buffer.data(fpg0 + 404);
    const auto *fpg0_405 = buffer.data(fpg0 + 405);
    const auto *fpg0_406 = buffer.data(fpg0 + 406);
    const auto *fpg0_407 = buffer.data(fpg0 + 407);
    const auto *fpg0_408 = buffer.data(fpg0 + 408);
    const auto *fpg0_409 = buffer.data(fpg0 + 409);
    const auto *fpg0_410 = buffer.data(fpg0 + 410);
    const auto *fpg0_419 = buffer.data(fpg0 + 419);

    const auto *fpf_181 = buffer.data(fpf + 181);
    const auto *fpf_186 = buffer.data(fpf + 186);
    const auto *fpf_196 = buffer.data(fpf + 196);
    const auto *fpf_197 = buffer.data(fpf + 197);
    const auto *fpf_199 = buffer.data(fpf + 199);
    const auto *fpf_206 = buffer.data(fpf + 206);
    const auto *fpf_209 = buffer.data(fpf + 209);
    const auto *fpf_216 = buffer.data(fpf + 216);
    const auto *fpf_219 = buffer.data(fpf + 219);
    const auto *fpf_226 = buffer.data(fpf + 226);
    const auto *fpf_229 = buffer.data(fpf + 229);
    const auto *fpf_236 = buffer.data(fpf + 236);
    const auto *fpf_238 = buffer.data(fpf + 238);
    const auto *fpf_239 = buffer.data(fpf + 239);
    const auto *fpf_246 = buffer.data(fpf + 246);
    const auto *fpf_249 = buffer.data(fpf + 249);
    const auto *fpf_258 = buffer.data(fpf + 258);
    const auto *fpf_259 = buffer.data(fpf + 259);
    const auto *fpf_266 = buffer.data(fpf + 266);
    const auto *fpf_268 = buffer.data(fpf + 268);
    const auto *fpf_269 = buffer.data(fpf + 269);
    const auto *fpf_270 = buffer.data(fpf + 270);
    const auto *fpf_271 = buffer.data(fpf + 271);
    const auto *fpf_272 = buffer.data(fpf + 272);
    const auto *fpf_279 = buffer.data(fpf + 279);

    const auto *fpg1_270 = buffer.data(fpg1 + 270);
    const auto *fpg1_271 = buffer.data(fpg1 + 271);
    const auto *fpg1_273 = buffer.data(fpg1 + 273);
    const auto *fpg1_274 = buffer.data(fpg1 + 274);
    const auto *fpg1_280 = buffer.data(fpg1 + 280);
    const auto *fpg1_285 = buffer.data(fpg1 + 285);
    const auto *fpg1_286 = buffer.data(fpg1 + 286);
    const auto *fpg1_288 = buffer.data(fpg1 + 288);
    const auto *fpg1_295 = buffer.data(fpg1 + 295);
    const auto *fpg1_297 = buffer.data(fpg1 + 297);
    const auto *fpg1_340 = buffer.data(fpg1 + 340);
    const auto *fpg1_359 = buffer.data(fpg1 + 359);
    const auto *fpg1_404 = buffer.data(fpg1 + 404);
    const auto *fpg1_405 = buffer.data(fpg1 + 405);
    const auto *fpg1_406 = buffer.data(fpg1 + 406);
    const auto *fpg1_407 = buffer.data(fpg1 + 407);
    const auto *fpg1_408 = buffer.data(fpg1 + 408);
    const auto *fpg1_409 = buffer.data(fpg1 + 409);
    const auto *fpg1_410 = buffer.data(fpg1 + 410);
    const auto *fpg1_419 = buffer.data(fpg1 + 419);

    const auto *gsg0_162 = buffer.data(gsg0 + 162);
    const auto *gsg0_164 = buffer.data(gsg0 + 164);
    const auto *gsg0_167 = buffer.data(gsg0 + 167);
    const auto *gsg0_170 = buffer.data(gsg0 + 170);
    const auto *gsg0_177 = buffer.data(gsg0 + 177);
    const auto *gsg0_179 = buffer.data(gsg0 + 179);
    const auto *gsg0_180 = buffer.data(gsg0 + 180);
    const auto *gsg0_181 = buffer.data(gsg0 + 181);
    const auto *gsg0_182 = buffer.data(gsg0 + 182);
    const auto *gsg0_183 = buffer.data(gsg0 + 183);
    const auto *gsg0_184 = buffer.data(gsg0 + 184);
    const auto *gsg0_185 = buffer.data(gsg0 + 185);
    const auto *gsg0_190 = buffer.data(gsg0 + 190);
    const auto *gsg0_192 = buffer.data(gsg0 + 192);
    const auto *gsg0_194 = buffer.data(gsg0 + 194);
    const auto *gsg0_205 = buffer.data(gsg0 + 205);
    const auto *gsg0_207 = buffer.data(gsg0 + 207);

    const auto *gsf_107 = buffer.data(gsf + 107);
    const auto *gsf_109 = buffer.data(gsf + 109);
    const auto *gsf_112 = buffer.data(gsf + 112);
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
    const auto *gsf_136 = buffer.data(gsf + 136);
    const auto *gsf_137 = buffer.data(gsf + 137);
    const auto *gsf_138 = buffer.data(gsf + 138);
    const auto *gsf_139 = buffer.data(gsf + 139);

    const auto *gsg1_162 = buffer.data(gsg1 + 162);
    const auto *gsg1_164 = buffer.data(gsg1 + 164);
    const auto *gsg1_167 = buffer.data(gsg1 + 167);
    const auto *gsg1_170 = buffer.data(gsg1 + 170);
    const auto *gsg1_177 = buffer.data(gsg1 + 177);
    const auto *gsg1_179 = buffer.data(gsg1 + 179);
    const auto *gsg1_180 = buffer.data(gsg1 + 180);
    const auto *gsg1_181 = buffer.data(gsg1 + 181);
    const auto *gsg1_182 = buffer.data(gsg1 + 182);
    const auto *gsg1_183 = buffer.data(gsg1 + 183);
    const auto *gsg1_184 = buffer.data(gsg1 + 184);
    const auto *gsg1_185 = buffer.data(gsg1 + 185);
    const auto *gsg1_190 = buffer.data(gsg1 + 190);
    const auto *gsg1_192 = buffer.data(gsg1 + 192);
    const auto *gsg1_194 = buffer.data(gsg1 + 194);
    const auto *gsg1_205 = buffer.data(gsg1 + 205);
    const auto *gsg1_207 = buffer.data(gsg1 + 207);

    const auto *gpd0_206 = buffer.data(gpd0 + 206);
    const auto *gpd0_208 = buffer.data(gpd0 + 208);
    const auto *gpd0_209 = buffer.data(gpd0 + 209);
    const auto *gpd0_210 = buffer.data(gpd0 + 210);
    const auto *gpd0_211 = buffer.data(gpd0 + 211);
    const auto *gpd0_212 = buffer.data(gpd0 + 212);
    const auto *gpd0_213 = buffer.data(gpd0 + 213);
    const auto *gpd0_214 = buffer.data(gpd0 + 214);
    const auto *gpd0_215 = buffer.data(gpd0 + 215);
    const auto *gpd0_222 = buffer.data(gpd0 + 222);
    const auto *gpd0_223 = buffer.data(gpd0 + 223);
    const auto *gpd0_224 = buffer.data(gpd0 + 224);
    const auto *gpd0_225 = buffer.data(gpd0 + 225);
    const auto *gpd0_226 = buffer.data(gpd0 + 226);
    const auto *gpd0_227 = buffer.data(gpd0 + 227);
    const auto *gpd0_228 = buffer.data(gpd0 + 228);
    const auto *gpd0_229 = buffer.data(gpd0 + 229);
    const auto *gpd0_230 = buffer.data(gpd0 + 230);
    const auto *gpd0_231 = buffer.data(gpd0 + 231);
    const auto *gpd0_232 = buffer.data(gpd0 + 232);
    const auto *gpd0_233 = buffer.data(gpd0 + 233);
    const auto *gpd0_240 = buffer.data(gpd0 + 240);
    const auto *gpd0_241 = buffer.data(gpd0 + 241);
    const auto *gpd0_242 = buffer.data(gpd0 + 242);
    const auto *gpd0_243 = buffer.data(gpd0 + 243);
    const auto *gpd0_244 = buffer.data(gpd0 + 244);
    const auto *gpd0_245 = buffer.data(gpd0 + 245);

    const auto *gpd1_206 = buffer.data(gpd1 + 206);
    const auto *gpd1_208 = buffer.data(gpd1 + 208);
    const auto *gpd1_209 = buffer.data(gpd1 + 209);
    const auto *gpd1_210 = buffer.data(gpd1 + 210);
    const auto *gpd1_211 = buffer.data(gpd1 + 211);
    const auto *gpd1_212 = buffer.data(gpd1 + 212);
    const auto *gpd1_213 = buffer.data(gpd1 + 213);
    const auto *gpd1_214 = buffer.data(gpd1 + 214);
    const auto *gpd1_215 = buffer.data(gpd1 + 215);
    const auto *gpd1_222 = buffer.data(gpd1 + 222);
    const auto *gpd1_223 = buffer.data(gpd1 + 223);
    const auto *gpd1_224 = buffer.data(gpd1 + 224);
    const auto *gpd1_225 = buffer.data(gpd1 + 225);
    const auto *gpd1_226 = buffer.data(gpd1 + 226);
    const auto *gpd1_227 = buffer.data(gpd1 + 227);
    const auto *gpd1_228 = buffer.data(gpd1 + 228);
    const auto *gpd1_229 = buffer.data(gpd1 + 229);
    const auto *gpd1_230 = buffer.data(gpd1 + 230);
    const auto *gpd1_231 = buffer.data(gpd1 + 231);
    const auto *gpd1_232 = buffer.data(gpd1 + 232);
    const auto *gpd1_233 = buffer.data(gpd1 + 233);
    const auto *gpd1_240 = buffer.data(gpd1 + 240);
    const auto *gpd1_241 = buffer.data(gpd1 + 241);
    const auto *gpd1_242 = buffer.data(gpd1 + 242);
    const auto *gpd1_243 = buffer.data(gpd1 + 243);
    const auto *gpd1_244 = buffer.data(gpd1 + 244);
    const auto *gpd1_245 = buffer.data(gpd1 + 245);

    const auto *gpf_329 = buffer.data(gpf + 329);
    const auto *gpf_336 = buffer.data(gpf + 336);
    const auto *gpf_337 = buffer.data(gpf + 337);
    const auto *gpf_338 = buffer.data(gpf + 338);
    const auto *gpf_339 = buffer.data(gpf + 339);
    const auto *gpf_342 = buffer.data(gpf + 342);
    const auto *gpf_344 = buffer.data(gpf + 344);
    const auto *gpf_345 = buffer.data(gpf + 345);
    const auto *gpf_346 = buffer.data(gpf + 346);
    const auto *gpf_347 = buffer.data(gpf + 347);
    const auto *gpf_348 = buffer.data(gpf + 348);
    const auto *gpf_349 = buffer.data(gpf + 349);
    const auto *gpf_350 = buffer.data(gpf + 350);
    const auto *gpf_351 = buffer.data(gpf + 351);
    const auto *gpf_352 = buffer.data(gpf + 352);
    const auto *gpf_353 = buffer.data(gpf + 353);
    const auto *gpf_354 = buffer.data(gpf + 354);
    const auto *gpf_355 = buffer.data(gpf + 355);
    const auto *gpf_356 = buffer.data(gpf + 356);
    const auto *gpf_357 = buffer.data(gpf + 357);
    const auto *gpf_358 = buffer.data(gpf + 358);
    const auto *gpf_359 = buffer.data(gpf + 359);
    const auto *gpf_366 = buffer.data(gpf + 366);
    const auto *gpf_367 = buffer.data(gpf + 367);
    const auto *gpf_368 = buffer.data(gpf + 368);
    const auto *gpf_369 = buffer.data(gpf + 369);
    const auto *gpf_370 = buffer.data(gpf + 370);
    const auto *gpf_371 = buffer.data(gpf + 371);
    const auto *gpf_372 = buffer.data(gpf + 372);
    const auto *gpf_373 = buffer.data(gpf + 373);
    const auto *gpf_374 = buffer.data(gpf + 374);
    const auto *gpf_375 = buffer.data(gpf + 375);
    const auto *gpf_376 = buffer.data(gpf + 376);
    const auto *gpf_377 = buffer.data(gpf + 377);
    const auto *gpf_378 = buffer.data(gpf + 378);
    const auto *gpf_379 = buffer.data(gpf + 379);
    const auto *gpf_380 = buffer.data(gpf + 380);
    const auto *gpf_381 = buffer.data(gpf + 381);
    const auto *gpf_382 = buffer.data(gpf + 382);
    const auto *gpf_383 = buffer.data(gpf + 383);
    const auto *gpf_384 = buffer.data(gpf + 384);
    const auto *gpf_385 = buffer.data(gpf + 385);
    const auto *gpf_386 = buffer.data(gpf + 386);
    const auto *gpf_387 = buffer.data(gpf + 387);
    const auto *gpf_388 = buffer.data(gpf + 388);
    const auto *gpf_389 = buffer.data(gpf + 389);
    const auto *gpf_396 = buffer.data(gpf + 396);
    const auto *gpf_397 = buffer.data(gpf + 397);
    const auto *gpf_398 = buffer.data(gpf + 398);
    const auto *gpf_399 = buffer.data(gpf + 399);
    const auto *gpf_400 = buffer.data(gpf + 400);
    const auto *gpf_401 = buffer.data(gpf + 401);
    const auto *gpf_402 = buffer.data(gpf + 402);
    const auto *gpf_403 = buffer.data(gpf + 403);
    const auto *gpf_404 = buffer.data(gpf + 404);
    const auto *gpf_405 = buffer.data(gpf + 405);
    const auto *gpf_406 = buffer.data(gpf + 406);
    const auto *gpf_407 = buffer.data(gpf + 407);
    const auto *gpf_408 = buffer.data(gpf + 408);
    const auto *gpf_409 = buffer.data(gpf + 409);

#pragma omp simd aligned(t_492, t_493, t_494, pb_z, pc_y, pc_z, fpf_209, gsg0_162, gsg0_164, \
                         gsf_107, gsf_109, gsg1_162, gsg1_164, \
                         gpf_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = pb_z[k] * gsg0_162[k]
                   + f_8 * gsf_107[k]
                   - f_7 * pc_z[k] * gsg1_162[k];

        t_493[k] = f_0 * fpf_209[k]
                   + f_4 * pc_y[k] * gpf_329[k];

        t_494[k] = pb_z[k] * gsg0_164[k]
                   + f_0 * gsf_109[k]
                   - f_7 * pc_z[k] * gsg1_164[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, pa_z, pb_x, pc_x, pc_z, fpg0_270, fpg0_271, \
                         fpg1_270, fpg1_271, gsg0_167, gsf_112, \
                         gsg1_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = pa_z[k] * fpg0_270[k]
                   - f_7 * pc_z[k] * fpg1_270[k];

        t_496[k] = pa_z[k] * fpg0_271[k]
                   - f_7 * pc_z[k] * fpg1_271[k];

        t_497[k] = pb_x[k] * gsg0_167[k]
                   + f_9 * gsf_112[k]
                   - f_7 * pc_x[k] * gsg1_167[k];
    }

#pragma omp simd aligned(t_498, t_499, t_500, pa_z, pb_x, pc_x, pc_z, fpg0_273, fpg0_274, \
                         fpf_181, fpg1_273, fpg1_274, gsg0_170, gsf_115, \
                         gsg1_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = pa_z[k] * fpg0_273[k]
                   - f_7 * pc_z[k] * fpg1_273[k];

        t_499[k] = pa_z[k] * fpg0_274[k]
                   + f_1 * fpf_181[k]
                   - f_7 * pc_z[k] * fpg1_274[k];

        t_500[k] = pb_x[k] * gsg0_170[k]
                   + f_8 * gsf_115[k]
                   - f_7 * pc_x[k] * gsg1_170[k];
    }

#pragma omp simd aligned(t_501, t_502, t_503, t_504, pc_x, gsf_116, gsf_117, gsf_118, gsf_119, \
                         gpf_336, gpf_337, gpf_338, gpf_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_501[k] = f_1 * gsf_116[k]
                   + f_4 * pc_x[k] * gpf_336[k];

        t_502[k] = f_1 * gsf_117[k]
                   + f_4 * pc_x[k] * gpf_337[k];

        t_503[k] = f_1 * gsf_118[k]
                   + f_4 * pc_x[k] * gpf_338[k];

        t_504[k] = f_1 * gsf_119[k]
                   + f_4 * pc_x[k] * gpf_339[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, pa_z, pb_x, pc_x, pc_z, fpg0_280, fpf_186, \
                         fpg1_280, gsg0_177, gsg1_177, gpf_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = pa_z[k] * fpg0_280[k]
                   - f_7 * pc_z[k] * fpg1_280[k];

        t_506[k] = f_1 * fpf_186[k]
                   + f_4 * pc_z[k] * gpf_336[k];

        t_507[k] = pb_x[k] * gsg0_177[k]
                   - f_7 * pc_x[k] * gsg1_177[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, pa_z, pb_x, pc_x, pc_y, pc_z, fpg0_285, fpf_219, \
                         fpg1_285, gsg0_179, gsg1_179, gpf_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = f_9 * fpf_219[k]
                   + f_4 * pc_y[k] * gpf_339[k];

        t_509[k] = pb_x[k] * gsg0_179[k]
                   - f_7 * pc_x[k] * gsg1_179[k];

        t_510[k] = pa_z[k] * fpg0_285[k]
                   - f_7 * pc_z[k] * fpg1_285[k];
    }

#pragma omp simd aligned(t_511, t_512, t_513, pa_z, pc_x, pc_z, fpg0_286, fpg0_288, fpg1_286, \
                         fpg1_288, gpd0_206, gpd1_206, gpf_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_511[k] = pa_z[k] * fpg0_286[k]
                   - f_7 * pc_z[k] * fpg1_286[k];

        t_512[k] = f_12 * gpd0_206[k]
                   - f_13 * gpd1_206[k]
                   + f_4 * pc_x[k] * gpf_342[k];

        t_513[k] = pa_z[k] * fpg0_288[k]
                   - f_7 * pc_z[k] * fpg1_288[k];
    }

#pragma omp simd aligned(t_514, t_515, t_516, t_517, t_518, pc_x, gpd0_208, gpd0_209, \
                         gpd1_208, gpd1_209, gpf_344, gpf_345, gpf_346, gpf_347, \
                         gpf_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_514[k] = f_5 * gpd0_208[k]
                   - f_6 * gpd1_208[k]
                   + f_4 * pc_x[k] * gpf_344[k];

        t_515[k] = f_5 * gpd0_209[k]
                   - f_6 * gpd1_209[k]
                   + f_4 * pc_x[k] * gpf_345[k];

        t_516[k] = f_4 * pc_x[k] * gpf_346[k];

        t_517[k] = f_4 * pc_x[k] * gpf_347[k];

        t_518[k] = f_4 * pc_x[k] * gpf_348[k];
    }

#pragma omp simd aligned(t_519, t_520, t_521, t_522, pa_z, pc_x, pc_z, fpg0_295, fpg0_297, \
                         fpf_196, fpf_197, fpg1_295, fpg1_297, gpf_346, \
                         gpf_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_519[k] = f_4 * pc_x[k] * gpf_349[k];

        t_520[k] = pa_z[k] * fpg0_295[k]
                   - f_7 * pc_z[k] * fpg1_295[k];

        t_521[k] = f_1 * fpf_196[k]
                   + f_4 * pc_z[k] * gpf_346[k];

        t_522[k] = pa_z[k] * fpg0_297[k]
                   + f_8 * fpf_197[k]
                   - f_7 * pc_z[k] * fpg1_297[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, pc_x, pc_y, pc_z, fpf_199, fpf_229, gsf_119, \
                         gpd0_209, gpd0_210, gpd1_209, gpd1_210, gpf_349, \
                         gpf_350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = f_9 * fpf_229[k]
                   + f_1 * gsf_119[k]
                   + f_4 * pc_y[k] * gpf_349[k];

        t_524[k] = f_1 * fpf_199[k]
                   + f_2 * gpd0_209[k]
                   - f_3 * gpd1_209[k]
                   + f_4 * pc_z[k] * gpf_349[k];

        t_525[k] = f_2 * gpd0_210[k]
                   - f_3 * gpd1_210[k]
                   + f_4 * pc_x[k] * gpf_350[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, pc_x, gpd0_211, gpd0_212, gpd0_213, gpd1_211, \
                         gpd1_212, gpd1_213, gpf_351, gpf_352, \
                         gpf_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_12 * gpd0_211[k]
                   - f_13 * gpd1_211[k]
                   + f_4 * pc_x[k] * gpf_351[k];

        t_527[k] = f_12 * gpd0_212[k]
                   - f_13 * gpd1_212[k]
                   + f_4 * pc_x[k] * gpf_352[k];

        t_528[k] = f_5 * gpd0_213[k]
                   - f_6 * gpd1_213[k]
                   + f_4 * pc_x[k] * gpf_353[k];
    }

#pragma omp simd aligned(t_529, t_530, t_531, t_532, t_533, pc_x, gpd0_214, gpd0_215, \
                         gpd1_214, gpd1_215, gpf_354, gpf_355, gpf_356, gpf_357, \
                         gpf_358 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_529[k] = f_5 * gpd0_214[k]
                   - f_6 * gpd1_214[k]
                   + f_4 * pc_x[k] * gpf_354[k];

        t_530[k] = f_5 * gpd0_215[k]
                   - f_6 * gpd1_215[k]
                   + f_4 * pc_x[k] * gpf_355[k];

        t_531[k] = f_4 * pc_x[k] * gpf_356[k];

        t_532[k] = f_4 * pc_x[k] * gpf_357[k];

        t_533[k] = f_4 * pc_x[k] * gpf_358[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, pc_x, pc_y, pc_z, fpf_206, fpf_236, gsf_116, \
                         gpd0_213, gpd1_213, gpf_356, gpf_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = f_4 * pc_x[k] * gpf_359[k];

        t_535[k] = f_9 * fpf_236[k]
                   + f_2 * gpd0_213[k]
                   - f_3 * gpd1_213[k]
                   + f_4 * pc_y[k] * gpf_356[k];

        t_536[k] = f_1 * fpf_206[k]
                   + f_1 * gsf_116[k]
                   + f_4 * pc_z[k] * gpf_356[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, pa_y, pc_y, dpg0_224, dpg1_224, fpg0_359, \
                         fpf_238, fpf_239, fpg1_359, gpd0_215, gpd1_215, gpf_358, \
                         gpf_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = f_9 * fpf_238[k]
                   + f_5 * gpd0_215[k]
                   - f_6 * gpd1_215[k]
                   + f_4 * pc_y[k] * gpf_358[k];

        t_538[k] = f_9 * fpf_239[k]
                   + f_4 * pc_y[k] * gpf_359[k];

        t_539[k] = f_10 * dpg0_224[k]
                   - f_11 * dpg1_224[k]
                   + pa_y[k] * fpg0_359[k]
                   - f_7 * pc_y[k] * fpg1_359[k];
    }

#pragma omp simd aligned(t_540, t_541, t_542, pb_x, pc_x, gsg0_180, gsg0_181, gsg0_182, \
                         gsf_120, gsf_121, gsf_122, gsg1_180, gsg1_181, \
                         gsg1_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_540[k] = pb_x[k] * gsg0_180[k]
                   + f_0 * gsf_120[k]
                   - f_7 * pc_x[k] * gsg1_180[k];

        t_541[k] = pb_x[k] * gsg0_181[k]
                   + f_9 * gsf_121[k]
                   - f_7 * pc_x[k] * gsg1_181[k];

        t_542[k] = pb_x[k] * gsg0_182[k]
                   + f_9 * gsf_122[k]
                   - f_7 * pc_x[k] * gsg1_182[k];
    }

#pragma omp simd aligned(t_543, t_544, t_545, pb_x, pc_x, gsg0_183, gsg0_184, gsg0_185, \
                         gsf_123, gsf_124, gsf_125, gsg1_183, gsg1_184, \
                         gsg1_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_543[k] = pb_x[k] * gsg0_183[k]
                   + f_8 * gsf_123[k]
                   - f_7 * pc_x[k] * gsg1_183[k];

        t_544[k] = pb_x[k] * gsg0_184[k]
                   + f_8 * gsf_124[k]
                   - f_7 * pc_x[k] * gsg1_184[k];

        t_545[k] = pb_x[k] * gsg0_185[k]
                   + f_8 * gsf_125[k]
                   - f_7 * pc_x[k] * gsg1_185[k];
    }

#pragma omp simd aligned(t_546, t_547, t_548, t_549, pc_x, gsf_126, gsf_127, gsf_128, gsf_129, \
                         gpf_366, gpf_367, gpf_368, gpf_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_546[k] = f_1 * gsf_126[k]
                   + f_4 * pc_x[k] * gpf_366[k];

        t_547[k] = f_1 * gsf_127[k]
                   + f_4 * pc_x[k] * gpf_367[k];

        t_548[k] = f_1 * gsf_128[k]
                   + f_4 * pc_x[k] * gpf_368[k];

        t_549[k] = f_1 * gsf_129[k]
                   + f_4 * pc_x[k] * gpf_369[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, t_553, pb_x, pc_x, pc_y, pc_z, fpf_216, fpf_249, \
                         gsg0_190, gsg0_192, gsg1_190, gsg1_192, gpf_366, \
                         gpf_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = pb_x[k] * gsg0_190[k]
                   - f_7 * pc_x[k] * gsg1_190[k];

        t_551[k] = f_8 * fpf_216[k]
                   + f_4 * pc_z[k] * gpf_366[k];

        t_552[k] = pb_x[k] * gsg0_192[k]
                   - f_7 * pc_x[k] * gsg1_192[k];

        t_553[k] = f_8 * fpf_249[k]
                   + f_4 * pc_y[k] * gpf_369[k];
    }

#pragma omp simd aligned(t_554, t_555, t_556, pb_x, pc_x, gsg0_194, gsg1_194, gpd0_222, \
                         gpd0_223, gpd1_222, gpd1_223, gpf_370, \
                         gpf_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_554[k] = pb_x[k] * gsg0_194[k]
                   - f_7 * pc_x[k] * gsg1_194[k];

        t_555[k] = f_2 * gpd0_222[k]
                   - f_3 * gpd1_222[k]
                   + f_4 * pc_x[k] * gpf_370[k];

        t_556[k] = f_12 * gpd0_223[k]
                   - f_13 * gpd1_223[k]
                   + f_4 * pc_x[k] * gpf_371[k];
    }

#pragma omp simd aligned(t_557, t_558, t_559, pc_x, gpd0_224, gpd0_225, gpd0_226, gpd1_224, \
                         gpd1_225, gpd1_226, gpf_372, gpf_373, \
                         gpf_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_557[k] = f_12 * gpd0_224[k]
                   - f_13 * gpd1_224[k]
                   + f_4 * pc_x[k] * gpf_372[k];

        t_558[k] = f_5 * gpd0_225[k]
                   - f_6 * gpd1_225[k]
                   + f_4 * pc_x[k] * gpf_373[k];

        t_559[k] = f_5 * gpd0_226[k]
                   - f_6 * gpd1_226[k]
                   + f_4 * pc_x[k] * gpf_374[k];
    }

#pragma omp simd aligned(t_560, t_561, t_562, t_563, t_564, pc_x, gpd0_227, gpd1_227, gpf_375, \
                         gpf_376, gpf_377, gpf_378, gpf_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_560[k] = f_5 * gpd0_227[k]
                   - f_6 * gpd1_227[k]
                   + f_4 * pc_x[k] * gpf_375[k];

        t_561[k] = f_4 * pc_x[k] * gpf_376[k];

        t_562[k] = f_4 * pc_x[k] * gpf_377[k];

        t_563[k] = f_4 * pc_x[k] * gpf_378[k];

        t_564[k] = f_4 * pc_x[k] * gpf_379[k];
    }

#pragma omp simd aligned(t_565, t_566, pa_z, pc_z, dpg0_160, dpg1_160, fpg0_340, fpf_226, \
                         fpg1_340, gpf_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_565[k] = f_14 * dpg0_160[k]
                   - f_15 * dpg1_160[k]
                   + pa_z[k] * fpg0_340[k]
                   - f_7 * pc_z[k] * fpg1_340[k];

        t_566[k] = f_8 * fpf_226[k]
                   + f_4 * pc_z[k] * gpf_376[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, pc_y, pc_z, fpf_229, fpf_258, fpf_259, gsf_128, \
                         gsf_129, gpd0_227, gpd1_227, gpf_378, \
                         gpf_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = f_8 * fpf_258[k]
                   + f_1 * gsf_128[k]
                   + f_5 * gpd0_227[k]
                   - f_6 * gpd1_227[k]
                   + f_4 * pc_y[k] * gpf_378[k];

        t_568[k] = f_8 * fpf_259[k]
                   + f_1 * gsf_129[k]
                   + f_4 * pc_y[k] * gpf_379[k];

        t_569[k] = f_8 * fpf_229[k]
                   + f_2 * gpd0_227[k]
                   - f_3 * gpd1_227[k]
                   + f_4 * pc_z[k] * gpf_379[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, pc_x, gpd0_228, gpd0_229, gpd0_230, gpd1_228, \
                         gpd1_229, gpd1_230, gpf_380, gpf_381, \
                         gpf_382 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = f_2 * gpd0_228[k]
                   - f_3 * gpd1_228[k]
                   + f_4 * pc_x[k] * gpf_380[k];

        t_571[k] = f_12 * gpd0_229[k]
                   - f_13 * gpd1_229[k]
                   + f_4 * pc_x[k] * gpf_381[k];

        t_572[k] = f_12 * gpd0_230[k]
                   - f_13 * gpd1_230[k]
                   + f_4 * pc_x[k] * gpf_382[k];
    }

#pragma omp simd aligned(t_573, t_574, t_575, t_576, pc_x, gpd0_231, gpd0_232, gpd0_233, \
                         gpd1_231, gpd1_232, gpd1_233, gpf_383, gpf_384, gpf_385, \
                         gpf_386 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_573[k] = f_5 * gpd0_231[k]
                   - f_6 * gpd1_231[k]
                   + f_4 * pc_x[k] * gpf_383[k];

        t_574[k] = f_5 * gpd0_232[k]
                   - f_6 * gpd1_232[k]
                   + f_4 * pc_x[k] * gpf_384[k];

        t_575[k] = f_5 * gpd0_233[k]
                   - f_6 * gpd1_233[k]
                   + f_4 * pc_x[k] * gpf_385[k];

        t_576[k] = f_4 * pc_x[k] * gpf_386[k];
    }

#pragma omp simd aligned(t_577, t_578, t_579, t_580, pc_x, pc_y, fpf_266, gpd0_231, gpd1_231, \
                         gpf_386, gpf_387, gpf_388, gpf_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_577[k] = f_4 * pc_x[k] * gpf_387[k];

        t_578[k] = f_4 * pc_x[k] * gpf_388[k];

        t_579[k] = f_4 * pc_x[k] * gpf_389[k];

        t_580[k] = f_8 * fpf_266[k]
                   + f_2 * gpd0_231[k]
                   - f_3 * gpd1_231[k]
                   + f_4 * pc_y[k] * gpf_386[k];
    }

#pragma omp simd aligned(t_581, t_582, t_583, pc_y, pc_z, fpf_236, fpf_268, fpf_269, gsf_126, \
                         gpd0_233, gpd1_233, gpf_386, gpf_388, \
                         gpf_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_581[k] = f_8 * fpf_236[k]
                   + f_1 * gsf_126[k]
                   + f_4 * pc_z[k] * gpf_386[k];

        t_582[k] = f_8 * fpf_268[k]
                   + f_5 * gpd0_233[k]
                   - f_6 * gpd1_233[k]
                   + f_4 * pc_y[k] * gpf_388[k];

        t_583[k] = f_8 * fpf_269[k]
                   + f_4 * pc_y[k] * gpf_389[k];
    }

#pragma omp simd aligned(t_584, t_585, t_586, pa_y, pc_y, dpg0_269, dpg1_269, fpg0_404, \
                         fpg0_405, fpg0_406, fpf_270, fpg1_404, fpg1_405, \
                         fpg1_406 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_584[k] = f_14 * dpg0_269[k]
                   - f_15 * dpg1_269[k]
                   + pa_y[k] * fpg0_404[k]
                   - f_7 * pc_y[k] * fpg1_404[k];

        t_585[k] = pa_y[k] * fpg0_405[k]
                   - f_7 * pc_y[k] * fpg1_405[k];

        t_586[k] = pa_y[k] * fpg0_406[k]
                   + f_1 * fpf_270[k]
                   - f_7 * pc_y[k] * fpg1_406[k];
    }

#pragma omp simd aligned(t_587, t_588, t_589, t_590, pa_y, pc_y, fpg0_407, fpg0_408, fpg0_409, \
                         fpg0_410, fpf_271, fpf_272, fpg1_407, fpg1_408, fpg1_409, \
                         fpg1_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_587[k] = pa_y[k] * fpg0_407[k]
                   - f_7 * pc_y[k] * fpg1_407[k];

        t_588[k] = pa_y[k] * fpg0_408[k]
                   + f_8 * fpf_271[k]
                   - f_7 * pc_y[k] * fpg1_408[k];

        t_589[k] = pa_y[k] * fpg0_409[k]
                   + f_1 * fpf_272[k]
                   - f_7 * pc_y[k] * fpg1_409[k];

        t_590[k] = pa_y[k] * fpg0_410[k]
                   - f_7 * pc_y[k] * fpg1_410[k];
    }

#pragma omp simd aligned(t_591, t_592, t_593, t_594, pc_x, gsf_136, gsf_137, gsf_138, gsf_139, \
                         gpf_396, gpf_397, gpf_398, gpf_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_591[k] = f_1 * gsf_136[k]
                   + f_4 * pc_x[k] * gpf_396[k];

        t_592[k] = f_1 * gsf_137[k]
                   + f_4 * pc_x[k] * gpf_397[k];

        t_593[k] = f_1 * gsf_138[k]
                   + f_4 * pc_x[k] * gpf_398[k];

        t_594[k] = f_1 * gsf_139[k]
                   + f_4 * pc_x[k] * gpf_399[k];
    }

#pragma omp simd aligned(t_595, t_596, t_597, t_598, pb_x, pc_x, pc_y, pc_z, fpf_246, fpf_279, \
                         gsg0_205, gsg0_207, gsg1_205, gsg1_207, gpf_396, \
                         gpf_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_595[k] = pb_x[k] * gsg0_205[k]
                   - f_7 * pc_x[k] * gsg1_205[k];

        t_596[k] = f_9 * fpf_246[k]
                   + f_4 * pc_z[k] * gpf_396[k];

        t_597[k] = pb_x[k] * gsg0_207[k]
                   - f_7 * pc_x[k] * gsg1_207[k];

        t_598[k] = f_1 * fpf_279[k]
                   + f_4 * pc_y[k] * gpf_399[k];
    }

#pragma omp simd aligned(t_599, t_600, t_601, pa_y, pc_x, pc_y, fpg0_419, fpg1_419, gpd0_240, \
                         gpd0_241, gpd1_240, gpd1_241, gpf_400, \
                         gpf_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_599[k] = pa_y[k] * fpg0_419[k]
                   - f_7 * pc_y[k] * fpg1_419[k];

        t_600[k] = f_2 * gpd0_240[k]
                   - f_3 * gpd1_240[k]
                   + f_4 * pc_x[k] * gpf_400[k];

        t_601[k] = f_12 * gpd0_241[k]
                   - f_13 * gpd1_241[k]
                   + f_4 * pc_x[k] * gpf_401[k];
    }

#pragma omp simd aligned(t_602, t_603, t_604, pc_x, gpd0_242, gpd0_243, gpd0_244, gpd1_242, \
                         gpd1_243, gpd1_244, gpf_402, gpf_403, \
                         gpf_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_602[k] = f_12 * gpd0_242[k]
                   - f_13 * gpd1_242[k]
                   + f_4 * pc_x[k] * gpf_402[k];

        t_603[k] = f_5 * gpd0_243[k]
                   - f_6 * gpd1_243[k]
                   + f_4 * pc_x[k] * gpf_403[k];

        t_604[k] = f_5 * gpd0_244[k]
                   - f_6 * gpd1_244[k]
                   + f_4 * pc_x[k] * gpf_404[k];
    }

#pragma omp simd aligned(t_605, t_606, t_607, t_608, t_609, pc_x, gpd0_245, gpd1_245, gpf_405, \
                         gpf_406, gpf_407, gpf_408, gpf_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_605[k] = f_5 * gpd0_245[k]
                   - f_6 * gpd1_245[k]
                   + f_4 * pc_x[k] * gpf_405[k];

        t_606[k] = f_4 * pc_x[k] * gpf_406[k];

        t_607[k] = f_4 * pc_x[k] * gpf_407[k];

        t_608[k] = f_4 * pc_x[k] * gpf_408[k];

        t_609[k] = f_4 * pc_x[k] * gpf_409[k];
    }
}

static auto
compute_prim_gpg_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t fpg0, const size_t fpf,
                                                          const size_t fpg1, const size_t gsg0,
                                                          const size_t gsf, const size_t gsg1,
                                                          const size_t gpd0, const size_t gpd1,
                                                          const size_t gpf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 1.5 / gamma;
    const auto f_3 = 1.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = gamma / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 1.5 / q;
    const auto f_12 = 1.0 / gamma;
    const auto f_13 = p / (gamma * q);

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *fpg0_435 = buffer.data(fpg0 + 435);
    const auto *fpg0_437 = buffer.data(fpg0 + 437);
    const auto *fpg0_440 = buffer.data(fpg0 + 440);
    const auto *fpg0_445 = buffer.data(fpg0 + 445);
    const auto *fpg0_447 = buffer.data(fpg0 + 447);
    const auto *fpg0_449 = buffer.data(fpg0 + 449);

    const auto *fpf_256 = buffer.data(fpf + 256);
    const auto *fpf_259 = buffer.data(fpf + 259);
    const auto *fpf_266 = buffer.data(fpf + 266);
    const auto *fpf_286 = buffer.data(fpf + 286);
    const auto *fpf_288 = buffer.data(fpf + 288);
    const auto *fpf_289 = buffer.data(fpf + 289);
    const auto *fpf_296 = buffer.data(fpf + 296);
    const auto *fpf_298 = buffer.data(fpf + 298);
    const auto *fpf_299 = buffer.data(fpf + 299);

    const auto *fpg1_435 = buffer.data(fpg1 + 435);
    const auto *fpg1_437 = buffer.data(fpg1 + 437);
    const auto *fpg1_440 = buffer.data(fpg1 + 440);
    const auto *fpg1_445 = buffer.data(fpg1 + 445);
    const auto *fpg1_447 = buffer.data(fpg1 + 447);
    const auto *fpg1_449 = buffer.data(fpg1 + 449);

    const auto *gsg0_210 = buffer.data(gsg0 + 210);
    const auto *gsg0_212 = buffer.data(gsg0 + 212);
    const auto *gsg0_213 = buffer.data(gsg0 + 213);
    const auto *gsg0_215 = buffer.data(gsg0 + 215);
    const auto *gsg0_220 = buffer.data(gsg0 + 220);
    const auto *gsg0_221 = buffer.data(gsg0 + 221);
    const auto *gsg0_222 = buffer.data(gsg0 + 222);
    const auto *gsg0_224 = buffer.data(gsg0 + 224);

    const auto *gsf_136 = buffer.data(gsf + 136);
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

    const auto *gsg1_210 = buffer.data(gsg1 + 210);
    const auto *gsg1_212 = buffer.data(gsg1 + 212);
    const auto *gsg1_213 = buffer.data(gsg1 + 213);
    const auto *gsg1_215 = buffer.data(gsg1 + 215);
    const auto *gsg1_220 = buffer.data(gsg1 + 220);
    const auto *gsg1_221 = buffer.data(gsg1 + 221);
    const auto *gsg1_222 = buffer.data(gsg1 + 222);
    const auto *gsg1_224 = buffer.data(gsg1 + 224);

    const auto *gpd0_243 = buffer.data(gpd0 + 243);
    const auto *gpd0_245 = buffer.data(gpd0 + 245);
    const auto *gpd0_247 = buffer.data(gpd0 + 247);
    const auto *gpd0_249 = buffer.data(gpd0 + 249);
    const auto *gpd0_250 = buffer.data(gpd0 + 250);
    const auto *gpd0_261 = buffer.data(gpd0 + 261);
    const auto *gpd0_264 = buffer.data(gpd0 + 264);
    const auto *gpd0_266 = buffer.data(gpd0 + 266);
    const auto *gpd0_267 = buffer.data(gpd0 + 267);
    const auto *gpd0_268 = buffer.data(gpd0 + 268);
    const auto *gpd0_269 = buffer.data(gpd0 + 269);

    const auto *gpd1_243 = buffer.data(gpd1 + 243);
    const auto *gpd1_245 = buffer.data(gpd1 + 245);
    const auto *gpd1_247 = buffer.data(gpd1 + 247);
    const auto *gpd1_249 = buffer.data(gpd1 + 249);
    const auto *gpd1_250 = buffer.data(gpd1 + 250);
    const auto *gpd1_261 = buffer.data(gpd1 + 261);
    const auto *gpd1_264 = buffer.data(gpd1 + 264);
    const auto *gpd1_266 = buffer.data(gpd1 + 266);
    const auto *gpd1_267 = buffer.data(gpd1 + 267);
    const auto *gpd1_268 = buffer.data(gpd1 + 268);
    const auto *gpd1_269 = buffer.data(gpd1 + 269);

    const auto *gpf_406 = buffer.data(gpf + 406);
    const auto *gpf_408 = buffer.data(gpf + 408);
    const auto *gpf_409 = buffer.data(gpf + 409);
    const auto *gpf_411 = buffer.data(gpf + 411);
    const auto *gpf_413 = buffer.data(gpf + 413);
    const auto *gpf_414 = buffer.data(gpf + 414);
    const auto *gpf_416 = buffer.data(gpf + 416);
    const auto *gpf_417 = buffer.data(gpf + 417);
    const auto *gpf_418 = buffer.data(gpf + 418);
    const auto *gpf_419 = buffer.data(gpf + 419);
    const auto *gpf_420 = buffer.data(gpf + 420);
    const auto *gpf_422 = buffer.data(gpf + 422);
    const auto *gpf_426 = buffer.data(gpf + 426);
    const auto *gpf_427 = buffer.data(gpf + 427);
    const auto *gpf_428 = buffer.data(gpf + 428);
    const auto *gpf_429 = buffer.data(gpf + 429);
    const auto *gpf_430 = buffer.data(gpf + 430);
    const auto *gpf_432 = buffer.data(gpf + 432);
    const auto *gpf_433 = buffer.data(gpf + 433);
    const auto *gpf_436 = buffer.data(gpf + 436);
    const auto *gpf_437 = buffer.data(gpf + 437);
    const auto *gpf_438 = buffer.data(gpf + 438);
    const auto *gpf_439 = buffer.data(gpf + 439);
    const auto *gpf_440 = buffer.data(gpf + 440);
    const auto *gpf_442 = buffer.data(gpf + 442);
    const auto *gpf_443 = buffer.data(gpf + 443);
    const auto *gpf_445 = buffer.data(gpf + 445);
    const auto *gpf_446 = buffer.data(gpf + 446);
    const auto *gpf_447 = buffer.data(gpf + 447);
    const auto *gpf_448 = buffer.data(gpf + 448);
    const auto *gpf_449 = buffer.data(gpf + 449);

#pragma omp simd aligned(t_610, t_611, t_612, pc_y, pc_z, fpf_256, fpf_286, fpf_288, gsf_136, \
                         gsf_138, gpd0_243, gpd0_245, gpd1_243, gpd1_245, gpf_406, \
                         gpf_408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = f_1 * fpf_286[k]
                   + f_1 * gsf_136[k]
                   + f_2 * gpd0_243[k]
                   - f_3 * gpd1_243[k]
                   + f_4 * pc_y[k] * gpf_406[k];

        t_611[k] = f_9 * fpf_256[k]
                   + f_4 * pc_z[k] * gpf_406[k];

        t_612[k] = f_1 * fpf_288[k]
                   + f_1 * gsf_138[k]
                   + f_5 * gpd0_245[k]
                   - f_6 * gpd1_245[k]
                   + f_4 * pc_y[k] * gpf_408[k];
    }

#pragma omp simd aligned(t_613, t_614, t_615, pa_y, pc_y, pc_z, fpg0_435, fpf_259, fpf_289, \
                         fpg1_435, gsf_139, gpd0_245, gpd1_245, \
                         gpf_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_613[k] = f_1 * fpf_289[k]
                   + f_1 * gsf_139[k]
                   + f_4 * pc_y[k] * gpf_409[k];

        t_614[k] = f_9 * fpf_259[k]
                   + f_2 * gpd0_245[k]
                   - f_3 * gpd1_245[k]
                   + f_4 * pc_z[k] * gpf_409[k];

        t_615[k] = pa_y[k] * fpg0_435[k]
                   - f_7 * pc_y[k] * fpg1_435[k];
    }

#pragma omp simd aligned(t_616, t_617, t_618, pa_y, pc_x, pc_y, fpg0_437, fpg1_437, gpd0_247, \
                         gpd0_249, gpd1_247, gpd1_249, gpf_411, \
                         gpf_413 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_616[k] = f_12 * gpd0_247[k]
                   - f_13 * gpd1_247[k]
                   + f_4 * pc_x[k] * gpf_411[k];

        t_617[k] = pa_y[k] * fpg0_437[k]
                   - f_7 * pc_y[k] * fpg1_437[k];

        t_618[k] = f_5 * gpd0_249[k]
                   - f_6 * gpd1_249[k]
                   + f_4 * pc_x[k] * gpf_413[k];
    }

#pragma omp simd aligned(t_619, t_620, t_621, t_622, t_623, pa_y, pc_x, pc_y, fpg0_440, \
                         fpg1_440, gpd0_250, gpd1_250, gpf_414, gpf_416, gpf_417, \
                         gpf_418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_619[k] = f_5 * gpd0_250[k]
                   - f_6 * gpd1_250[k]
                   + f_4 * pc_x[k] * gpf_414[k];

        t_620[k] = pa_y[k] * fpg0_440[k]
                   - f_7 * pc_y[k] * fpg1_440[k];

        t_621[k] = f_4 * pc_x[k] * gpf_416[k];

        t_622[k] = f_4 * pc_x[k] * gpf_417[k];

        t_623[k] = f_4 * pc_x[k] * gpf_418[k];
    }

#pragma omp simd aligned(t_624, t_625, t_626, pa_y, pc_x, pc_y, pc_z, fpg0_445, fpf_266, \
                         fpf_296, fpg1_445, gsf_136, gpf_416, gpf_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_624[k] = f_4 * pc_x[k] * gpf_419[k];

        t_625[k] = pa_y[k] * fpg0_445[k]
                   + f_0 * fpf_296[k]
                   - f_7 * pc_y[k] * fpg1_445[k];

        t_626[k] = f_9 * fpf_266[k]
                   + f_1 * gsf_136[k]
                   + f_4 * pc_z[k] * gpf_416[k];
    }

#pragma omp simd aligned(t_627, t_628, t_629, pa_y, pc_y, fpg0_447, fpg0_449, fpf_298, \
                         fpf_299, fpg1_447, fpg1_449, gpf_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_627[k] = pa_y[k] * fpg0_447[k]
                   + f_8 * fpf_298[k]
                   - f_7 * pc_y[k] * fpg1_447[k];

        t_628[k] = f_1 * fpf_299[k]
                   + f_4 * pc_y[k] * gpf_419[k];

        t_629[k] = pa_y[k] * fpg0_449[k]
                   - f_7 * pc_y[k] * fpg1_449[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, pb_x, pc_x, pc_y, gsg0_210, gsg0_212, gsf_140, \
                         gsf_142, gsg1_210, gsg1_212, gpf_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = pb_x[k] * gsg0_210[k]
                   + f_0 * gsf_140[k]
                   - f_7 * pc_x[k] * gsg1_210[k];

        t_631[k] = f_4 * pc_y[k] * gpf_420[k];

        t_632[k] = pb_x[k] * gsg0_212[k]
                   + f_9 * gsf_142[k]
                   - f_7 * pc_x[k] * gsg1_212[k];
    }

#pragma omp simd aligned(t_633, t_634, t_635, t_636, pb_x, pc_x, pc_y, gsg0_213, gsg0_215, \
                         gsf_143, gsf_145, gsf_146, gsg1_213, gsg1_215, gpf_422, \
                         gpf_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_633[k] = pb_x[k] * gsg0_213[k]
                   + f_8 * gsf_143[k]
                   - f_7 * pc_x[k] * gsg1_213[k];

        t_634[k] = f_4 * pc_y[k] * gpf_422[k];

        t_635[k] = pb_x[k] * gsg0_215[k]
                   + f_8 * gsf_145[k]
                   - f_7 * pc_x[k] * gsg1_215[k];

        t_636[k] = f_1 * gsf_146[k]
                   + f_4 * pc_x[k] * gpf_426[k];
    }

#pragma omp simd aligned(t_637, t_638, t_639, t_640, pb_x, pc_x, gsg0_220, gsf_147, gsf_148, \
                         gsf_149, gsg1_220, gpf_427, gpf_428, gpf_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_637[k] = f_1 * gsf_147[k]
                   + f_4 * pc_x[k] * gpf_427[k];

        t_638[k] = f_1 * gsf_148[k]
                   + f_4 * pc_x[k] * gpf_428[k];

        t_639[k] = f_1 * gsf_149[k]
                   + f_4 * pc_x[k] * gpf_429[k];

        t_640[k] = pb_x[k] * gsg0_220[k]
                   - f_7 * pc_x[k] * gsg1_220[k];
    }

#pragma omp simd aligned(t_641, t_642, t_643, t_644, pb_x, pc_x, pc_y, gsg0_221, gsg0_222, \
                         gsg0_224, gsg1_221, gsg1_222, gsg1_224, \
                         gpf_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_641[k] = pb_x[k] * gsg0_221[k]
                   - f_7 * pc_x[k] * gsg1_221[k];

        t_642[k] = pb_x[k] * gsg0_222[k]
                   - f_7 * pc_x[k] * gsg1_222[k];

        t_643[k] = f_4 * pc_y[k] * gpf_429[k];

        t_644[k] = pb_x[k] * gsg0_224[k]
                   - f_7 * pc_x[k] * gsg1_224[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, t_648, pb_y, pc_x, pc_y, gsg0_210, gsg0_212, \
                         gsf_140, gsg1_210, gsg1_212, gpd0_261, gpd1_261, gpf_430, \
                         gpf_433 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = pb_y[k] * gsg0_210[k]
                   - f_7 * pc_y[k] * gsg1_210[k];

        t_646[k] = f_1 * gsf_140[k]
                   + f_4 * pc_y[k] * gpf_430[k];

        t_647[k] = pb_y[k] * gsg0_212[k]
                   - f_7 * pc_y[k] * gsg1_212[k];

        t_648[k] = f_5 * gpd0_261[k]
                   - f_6 * gpd1_261[k]
                   + f_4 * pc_x[k] * gpf_433[k];
    }

#pragma omp simd aligned(t_649, t_650, t_651, t_652, t_653, pb_y, pc_x, pc_y, gsg0_215, \
                         gsf_142, gsg1_215, gpf_432, gpf_436, gpf_437, \
                         gpf_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_649[k] = f_1 * gsf_142[k]
                   + f_4 * pc_y[k] * gpf_432[k];

        t_650[k] = pb_y[k] * gsg0_215[k]
                   - f_7 * pc_y[k] * gsg1_215[k];

        t_651[k] = f_4 * pc_x[k] * gpf_436[k];

        t_652[k] = f_4 * pc_x[k] * gpf_437[k];

        t_653[k] = f_4 * pc_x[k] * gpf_438[k];
    }

#pragma omp simd aligned(t_654, t_655, t_656, pb_y, pc_x, pc_y, gsg0_220, gsg0_221, gsf_146, \
                         gsf_147, gsg1_220, gsg1_221, gpf_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_654[k] = f_4 * pc_x[k] * gpf_439[k];

        t_655[k] = pb_y[k] * gsg0_220[k]
                   + f_0 * gsf_146[k]
                   - f_7 * pc_y[k] * gsg1_220[k];

        t_656[k] = pb_y[k] * gsg0_221[k]
                   + f_9 * gsf_147[k]
                   - f_7 * pc_y[k] * gsg1_221[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, pb_y, pc_y, gsg0_222, gsg0_224, gsf_148, \
                         gsf_149, gsg1_222, gsg1_224, gpf_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = pb_y[k] * gsg0_222[k]
                   + f_8 * gsf_148[k]
                   - f_7 * pc_y[k] * gsg1_222[k];

        t_658[k] = f_1 * gsf_149[k]
                   + f_4 * pc_y[k] * gpf_439[k];

        t_659[k] = pb_y[k] * gsg0_224[k]
                   - f_7 * pc_y[k] * gsg1_224[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, t_664, pc_x, pc_y, gpd0_264, gpd0_266, \
                         gpd0_267, gpd1_264, gpd1_266, gpd1_267, gpf_440, gpf_442, \
                         gpf_443 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = f_2 * gpd0_264[k]
                   - f_3 * gpd1_264[k]
                   + f_4 * pc_x[k] * gpf_440[k];

        t_661[k] = f_4 * pc_y[k] * gpf_440[k];

        t_662[k] = f_12 * gpd0_266[k]
                   - f_13 * gpd1_266[k]
                   + f_4 * pc_x[k] * gpf_442[k];

        t_663[k] = f_5 * gpd0_267[k]
                   - f_6 * gpd1_267[k]
                   + f_4 * pc_x[k] * gpf_443[k];

        t_664[k] = f_4 * pc_y[k] * gpf_442[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, t_668, t_669, pc_x, gpd0_269, gpd1_269, gpf_445, \
                         gpf_446, gpf_447, gpf_448, gpf_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = f_5 * gpd0_269[k]
                   - f_6 * gpd1_269[k]
                   + f_4 * pc_x[k] * gpf_445[k];

        t_666[k] = f_4 * pc_x[k] * gpf_446[k];

        t_667[k] = f_4 * pc_x[k] * gpf_447[k];

        t_668[k] = f_4 * pc_x[k] * gpf_448[k];

        t_669[k] = f_4 * pc_x[k] * gpf_449[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, t_673, pc_y, gpd0_267, gpd0_268, gpd0_269, \
                         gpd1_267, gpd1_268, gpd1_269, gpf_446, gpf_447, gpf_448, \
                         gpf_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_2 * gpd0_267[k]
                   - f_3 * gpd1_267[k]
                   + f_4 * pc_y[k] * gpf_446[k];

        t_671[k] = f_12 * gpd0_268[k]
                   - f_13 * gpd1_268[k]
                   + f_4 * pc_y[k] * gpf_447[k];

        t_672[k] = f_5 * gpd0_269[k]
                   - f_6 * gpd1_269[k]
                   + f_4 * pc_y[k] * gpf_448[k];

        t_673[k] = f_4 * pc_y[k] * gpf_449[k];
    }

#pragma omp simd aligned(t_674, pc_z, fpf_299, gsf_149, gpd0_269, gpd1_269, \
                         gpf_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_674[k] = f_0 * fpf_299[k]
                   + f_1 * gsf_149[k]
                   + f_2 * gpd0_269[k]
                   - f_3 * gpd1_269[k]
                   + f_4 * pc_z[k] * gpf_449[k];
    }
}

auto
compute_prim_gpg_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t dpg0,
                                                   const size_t dpg1, const size_t fpg0,
                                                   const size_t fpf, const size_t fpg1,
                                                   const size_t gsg0, const size_t gsf,
                                                   const size_t gsg1, const size_t gpd0,
                                                   const size_t gpd1, const size_t gpf,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_gpg_three_center_electron_repulsion_0_piece0(buffer, target, pa, pb, pc, dpg0,
                                                              dpg1, fpg0, fpf, fpg1, gsg0, gsf,
                                                              gsg1, gpd0, gpd1, gpf, ncols,
                                                              gamma, p, q);

    compute_prim_gpg_three_center_electron_repulsion_0_piece1(buffer, target, pa, pb, pc, dpg0,
                                                              dpg1, fpg0, fpf, fpg1, gsg0, gsf,
                                                              gsg1, gpd0, gpd1, gpf, ncols,
                                                              gamma, p, q);

    compute_prim_gpg_three_center_electron_repulsion_0_piece2(buffer, target, pa, pb, pc, dpg0,
                                                              dpg1, fpg0, fpf, fpg1, gsg0, gsf,
                                                              gsg1, gpd0, gpd1, gpf, ncols,
                                                              gamma, p, q);

    compute_prim_gpg_three_center_electron_repulsion_0_piece3(buffer, target, pa, pb, pc, fpg0,
                                                              fpf, fpg1, gsg0, gsf, gsg1, gpd0,
                                                              gpd1, gpf, ncols, gamma, p, q);

    compute_prim_gpg_three_center_electron_repulsion_0_piece4(buffer, target, pa, pb, pc, dpg0,
                                                              dpg1, fpg0, fpf, fpg1, gsg0, gsf,
                                                              gsg1, gpd0, gpd1, gpf, ncols,
                                                              gamma, p, q);

    compute_prim_gpg_three_center_electron_repulsion_0_piece5(buffer, target, pa, pb, pc, fpg0,
                                                              fpf, fpg1, gsg0, gsf, gsg1, gpd0,
                                                              gpd1, gpf, ncols, gamma, p, q);
}

}  // namespace simdt3ceri
