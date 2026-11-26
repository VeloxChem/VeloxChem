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


#include "SimdThreeCenterElectronRepulsionVrrRecGPH.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_gph_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t dph0, const size_t dph1,
                                                          const size_t fph0, const size_t fpg,
                                                          const size_t fph1, const size_t gsh0,
                                                          const size_t gsg, const size_t gsh1,
                                                          const size_t gpf0, const size_t gpf1,
                                                          const size_t gpg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 2.0 / gamma;
    const auto f_3 = 2.0 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = gamma / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 2.5 / q;
    const auto f_13 = 1.0 / p;
    const auto f_14 = gamma / (p * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dph0_99 = buffer.data(dph0 + 99);

    const auto *dph1_99 = buffer.data(dph1 + 99);

    const auto *fph0_0 = buffer.data(fph0 + 0);
    const auto *fph0_3 = buffer.data(fph0 + 3);
    const auto *fph0_5 = buffer.data(fph0 + 5);
    const auto *fph0_6 = buffer.data(fph0 + 6);
    const auto *fph0_9 = buffer.data(fph0 + 9);
    const auto *fph0_14 = buffer.data(fph0 + 14);
    const auto *fph0_20 = buffer.data(fph0 + 20);
    const auto *fph0_42 = buffer.data(fph0 + 42);
    const auto *fph0_47 = buffer.data(fph0 + 47);
    const auto *fph0_51 = buffer.data(fph0 + 51);
    const auto *fph0_62 = buffer.data(fph0 + 62);
    const auto *fph0_99 = buffer.data(fph0 + 99);

    const auto *fpg_0 = buffer.data(fpg + 0);
    const auto *fpg_1 = buffer.data(fpg + 1);
    const auto *fpg_3 = buffer.data(fpg + 3);
    const auto *fpg_5 = buffer.data(fpg + 5);
    const auto *fpg_10 = buffer.data(fpg + 10);
    const auto *fpg_12 = buffer.data(fpg + 12);
    const auto *fpg_14 = buffer.data(fpg + 14);
    const auto *fpg_15 = buffer.data(fpg + 15);
    const auto *fpg_20 = buffer.data(fpg + 20);
    const auto *fpg_25 = buffer.data(fpg + 25);
    const auto *fpg_27 = buffer.data(fpg + 27);
    const auto *fpg_29 = buffer.data(fpg + 29);
    const auto *fpg_30 = buffer.data(fpg + 30);
    const auto *fpg_35 = buffer.data(fpg + 35);
    const auto *fpg_40 = buffer.data(fpg + 40);
    const auto *fpg_42 = buffer.data(fpg + 42);
    const auto *fpg_44 = buffer.data(fpg + 44);
    const auto *fpg_55 = buffer.data(fpg + 55);
    const auto *fpg_57 = buffer.data(fpg + 57);
    const auto *fpg_58 = buffer.data(fpg + 58);
    const auto *fpg_60 = buffer.data(fpg + 60);
    const auto *fpg_63 = buffer.data(fpg + 63);
    const auto *fpg_66 = buffer.data(fpg + 66);
    const auto *fpg_70 = buffer.data(fpg + 70);
    const auto *fpg_72 = buffer.data(fpg + 72);
    const auto *fpg_73 = buffer.data(fpg + 73);
    const auto *fpg_74 = buffer.data(fpg + 74);
    const auto *fpg_85 = buffer.data(fpg + 85);
    const auto *fpg_87 = buffer.data(fpg + 87);
    const auto *fpg_88 = buffer.data(fpg + 88);
    const auto *fpg_89 = buffer.data(fpg + 89);

    const auto *fph1_0 = buffer.data(fph1 + 0);
    const auto *fph1_3 = buffer.data(fph1 + 3);
    const auto *fph1_5 = buffer.data(fph1 + 5);
    const auto *fph1_6 = buffer.data(fph1 + 6);
    const auto *fph1_9 = buffer.data(fph1 + 9);
    const auto *fph1_14 = buffer.data(fph1 + 14);
    const auto *fph1_20 = buffer.data(fph1 + 20);
    const auto *fph1_42 = buffer.data(fph1 + 42);
    const auto *fph1_47 = buffer.data(fph1 + 47);
    const auto *fph1_51 = buffer.data(fph1 + 51);
    const auto *fph1_62 = buffer.data(fph1 + 62);
    const auto *fph1_99 = buffer.data(fph1 + 99);

    const auto *gsh0_0 = buffer.data(gsh0 + 0);
    const auto *gsh0_3 = buffer.data(gsh0 + 3);
    const auto *gsh0_5 = buffer.data(gsh0 + 5);
    const auto *gsh0_6 = buffer.data(gsh0 + 6);
    const auto *gsh0_9 = buffer.data(gsh0 + 9);
    const auto *gsh0_15 = buffer.data(gsh0 + 15);
    const auto *gsh0_17 = buffer.data(gsh0 + 17);
    const auto *gsh0_18 = buffer.data(gsh0 + 18);
    const auto *gsh0_20 = buffer.data(gsh0 + 20);
    const auto *gsh0_24 = buffer.data(gsh0 + 24);
    const auto *gsh0_27 = buffer.data(gsh0 + 27);
    const auto *gsh0_36 = buffer.data(gsh0 + 36);
    const auto *gsh0_38 = buffer.data(gsh0 + 38);
    const auto *gsh0_39 = buffer.data(gsh0 + 39);

    const auto *gsg_0 = buffer.data(gsg + 0);
    const auto *gsg_1 = buffer.data(gsg + 1);
    const auto *gsg_2 = buffer.data(gsg + 2);
    const auto *gsg_3 = buffer.data(gsg + 3);
    const auto *gsg_5 = buffer.data(gsg + 5);
    const auto *gsg_6 = buffer.data(gsg + 6);
    const auto *gsg_9 = buffer.data(gsg + 9);
    const auto *gsg_10 = buffer.data(gsg + 10);
    const auto *gsg_12 = buffer.data(gsg + 12);
    const auto *gsg_13 = buffer.data(gsg + 13);
    const auto *gsg_14 = buffer.data(gsg + 14);
    const auto *gsg_15 = buffer.data(gsg + 15);
    const auto *gsg_16 = buffer.data(gsg + 16);
    const auto *gsg_18 = buffer.data(gsg + 18);
    const auto *gsg_20 = buffer.data(gsg + 20);
    const auto *gsg_21 = buffer.data(gsg + 21);
    const auto *gsg_25 = buffer.data(gsg + 25);
    const auto *gsg_26 = buffer.data(gsg + 26);
    const auto *gsg_27 = buffer.data(gsg + 27);
    const auto *gsg_28 = buffer.data(gsg + 28);
    const auto *gsg_29 = buffer.data(gsg + 29);

    const auto *gsh1_0 = buffer.data(gsh1 + 0);
    const auto *gsh1_3 = buffer.data(gsh1 + 3);
    const auto *gsh1_5 = buffer.data(gsh1 + 5);
    const auto *gsh1_6 = buffer.data(gsh1 + 6);
    const auto *gsh1_9 = buffer.data(gsh1 + 9);
    const auto *gsh1_15 = buffer.data(gsh1 + 15);
    const auto *gsh1_17 = buffer.data(gsh1 + 17);
    const auto *gsh1_18 = buffer.data(gsh1 + 18);
    const auto *gsh1_20 = buffer.data(gsh1 + 20);
    const auto *gsh1_24 = buffer.data(gsh1 + 24);
    const auto *gsh1_27 = buffer.data(gsh1 + 27);
    const auto *gsh1_36 = buffer.data(gsh1 + 36);
    const auto *gsh1_38 = buffer.data(gsh1 + 38);
    const auto *gsh1_39 = buffer.data(gsh1 + 39);

    const auto *gpf0_0 = buffer.data(gpf0 + 0);
    const auto *gpf0_1 = buffer.data(gpf0 + 1);
    const auto *gpf0_2 = buffer.data(gpf0 + 2);
    const auto *gpf0_6 = buffer.data(gpf0 + 6);
    const auto *gpf0_8 = buffer.data(gpf0 + 8);
    const auto *gpf0_9 = buffer.data(gpf0 + 9);
    const auto *gpf0_28 = buffer.data(gpf0 + 28);
    const auto *gpf0_29 = buffer.data(gpf0 + 29);
    const auto *gpf0_36 = buffer.data(gpf0 + 36);
    const auto *gpf0_37 = buffer.data(gpf0 + 37);
    const auto *gpf0_40 = buffer.data(gpf0 + 40);
    const auto *gpf0_42 = buffer.data(gpf0 + 42);
    const auto *gpf0_43 = buffer.data(gpf0 + 43);
    const auto *gpf0_46 = buffer.data(gpf0 + 46);
    const auto *gpf0_47 = buffer.data(gpf0 + 47);
    const auto *gpf0_49 = buffer.data(gpf0 + 49);

    const auto *gpf1_0 = buffer.data(gpf1 + 0);
    const auto *gpf1_1 = buffer.data(gpf1 + 1);
    const auto *gpf1_2 = buffer.data(gpf1 + 2);
    const auto *gpf1_6 = buffer.data(gpf1 + 6);
    const auto *gpf1_8 = buffer.data(gpf1 + 8);
    const auto *gpf1_9 = buffer.data(gpf1 + 9);
    const auto *gpf1_28 = buffer.data(gpf1 + 28);
    const auto *gpf1_29 = buffer.data(gpf1 + 29);
    const auto *gpf1_36 = buffer.data(gpf1 + 36);
    const auto *gpf1_37 = buffer.data(gpf1 + 37);
    const auto *gpf1_40 = buffer.data(gpf1 + 40);
    const auto *gpf1_42 = buffer.data(gpf1 + 42);
    const auto *gpf1_43 = buffer.data(gpf1 + 43);
    const auto *gpf1_46 = buffer.data(gpf1 + 46);
    const auto *gpf1_47 = buffer.data(gpf1 + 47);
    const auto *gpf1_49 = buffer.data(gpf1 + 49);

    const auto *gpg_0 = buffer.data(gpg + 0);
    const auto *gpg_1 = buffer.data(gpg + 1);
    const auto *gpg_2 = buffer.data(gpg + 2);
    const auto *gpg_3 = buffer.data(gpg + 3);
    const auto *gpg_5 = buffer.data(gpg + 5);
    const auto *gpg_6 = buffer.data(gpg + 6);
    const auto *gpg_9 = buffer.data(gpg + 9);
    const auto *gpg_10 = buffer.data(gpg + 10);
    const auto *gpg_12 = buffer.data(gpg + 12);
    const auto *gpg_13 = buffer.data(gpg + 13);
    const auto *gpg_14 = buffer.data(gpg + 14);
    const auto *gpg_15 = buffer.data(gpg + 15);
    const auto *gpg_17 = buffer.data(gpg + 17);
    const auto *gpg_18 = buffer.data(gpg + 18);
    const auto *gpg_20 = buffer.data(gpg + 20);
    const auto *gpg_21 = buffer.data(gpg + 21);
    const auto *gpg_24 = buffer.data(gpg + 24);
    const auto *gpg_25 = buffer.data(gpg + 25);
    const auto *gpg_27 = buffer.data(gpg + 27);
    const auto *gpg_29 = buffer.data(gpg + 29);
    const auto *gpg_30 = buffer.data(gpg + 30);
    const auto *gpg_32 = buffer.data(gpg + 32);
    const auto *gpg_33 = buffer.data(gpg + 33);
    const auto *gpg_35 = buffer.data(gpg + 35);
    const auto *gpg_36 = buffer.data(gpg + 36);
    const auto *gpg_39 = buffer.data(gpg + 39);
    const auto *gpg_40 = buffer.data(gpg + 40);
    const auto *gpg_42 = buffer.data(gpg + 42);
    const auto *gpg_43 = buffer.data(gpg + 43);
    const auto *gpg_44 = buffer.data(gpg + 44);
    const auto *gpg_45 = buffer.data(gpg + 45);
    const auto *gpg_46 = buffer.data(gpg + 46);
    const auto *gpg_48 = buffer.data(gpg + 48);
    const auto *gpg_50 = buffer.data(gpg + 50);
    const auto *gpg_51 = buffer.data(gpg + 51);
    const auto *gpg_55 = buffer.data(gpg + 55);
    const auto *gpg_56 = buffer.data(gpg + 56);
    const auto *gpg_57 = buffer.data(gpg + 57);
    const auto *gpg_58 = buffer.data(gpg + 58);
    const auto *gpg_59 = buffer.data(gpg + 59);
    const auto *gpg_60 = buffer.data(gpg + 60);
    const auto *gpg_61 = buffer.data(gpg + 61);
    const auto *gpg_62 = buffer.data(gpg + 62);
    const auto *gpg_63 = buffer.data(gpg + 63);
    const auto *gpg_65 = buffer.data(gpg + 65);
    const auto *gpg_66 = buffer.data(gpg + 66);
    const auto *gpg_70 = buffer.data(gpg + 70);
    const auto *gpg_71 = buffer.data(gpg + 71);
    const auto *gpg_72 = buffer.data(gpg + 72);
    const auto *gpg_73 = buffer.data(gpg + 73);
    const auto *gpg_74 = buffer.data(gpg + 74);
    const auto *gpg_75 = buffer.data(gpg + 75);
    const auto *gpg_76 = buffer.data(gpg + 76);
    const auto *gpg_78 = buffer.data(gpg + 78);
    const auto *gpg_80 = buffer.data(gpg + 80);
    const auto *gpg_81 = buffer.data(gpg + 81);
    const auto *gpg_85 = buffer.data(gpg + 85);
    const auto *gpg_87 = buffer.data(gpg + 87);
    const auto *gpg_88 = buffer.data(gpg + 88);
    const auto *gpg_89 = buffer.data(gpg + 89);
    const auto *gpg_90 = buffer.data(gpg + 90);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, fpg_0, gsg_0, gpf0_0, \
                         gpf1_0, gpg_0, gpg_1, gpg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fpg_0[k]
                 + f_1 * gsg_0[k]
                 + f_2 * gpf0_0[k]
                 - f_3 * gpf1_0[k]
                 + f_4 * pc_x[k] * gpg_0[k];

        t_1[k] = f_4 * pc_y[k] * gpg_0[k];

        t_2[k] = f_4 * pc_z[k] * gpg_0[k];

        t_3[k] = f_5 * gpf0_0[k]
                 - f_6 * gpf1_0[k]
                 + f_4 * pc_y[k] * gpg_1[k];

        t_4[k] = f_4 * pc_y[k] * gpg_2[k];

        t_5[k] = f_5 * gpf0_0[k]
                 - f_6 * gpf1_0[k]
                 + f_4 * pc_z[k] * gpg_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pc_y, pc_z, gpf0_1, gpf0_2, gpf1_1, gpf1_2, \
                         gpg_3, gpg_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_7 * gpf0_1[k]
                 - f_8 * gpf1_1[k]
                 + f_4 * pc_y[k] * gpg_3[k];

        t_7[k] = f_4 * pc_z[k] * gpg_3[k];

        t_8[k] = f_4 * pc_y[k] * gpg_5[k];

        t_9[k] = f_7 * gpf0_2[k]
                 - f_8 * gpf1_2[k]
                 + f_4 * pc_z[k] * gpg_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pc_x, pc_y, pc_z, fpg_10, fpg_12, gsg_10, \
                         gsg_12, gpg_6, gpg_9, gpg_10, gpg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * fpg_10[k]
                  + f_1 * gsg_10[k]
                  + f_4 * pc_x[k] * gpg_10[k];

        t_11[k] = f_4 * pc_z[k] * gpg_6[k];

        t_12[k] = f_0 * fpg_12[k]
                  + f_1 * gsg_12[k]
                  + f_4 * pc_x[k] * gpg_12[k];

        t_13[k] = f_4 * pc_y[k] * gpg_9[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pc_x, pc_y, pc_z, fpg_14, gsg_14, gpf0_6, \
                         gpf0_8, gpf1_6, gpf1_8, gpg_10, gpg_12, \
                         gpg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_0 * fpg_14[k]
                  + f_1 * gsg_14[k]
                  + f_4 * pc_x[k] * gpg_14[k];

        t_15[k] = f_2 * gpf0_6[k]
                  - f_3 * gpf1_6[k]
                  + f_4 * pc_y[k] * gpg_10[k];

        t_16[k] = f_4 * pc_z[k] * gpg_10[k];

        t_17[k] = f_7 * gpf0_8[k]
                  - f_8 * gpf1_8[k]
                  + f_4 * pc_y[k] * gpg_12[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pb_y, pc_y, pc_z, gsh0_0, gsg_0, \
                         gsh1_0, gpf0_9, gpf1_9, gpg_13, gpg_14, \
                         gpg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * gpf0_9[k]
                  - f_6 * gpf1_9[k]
                  + f_4 * pc_y[k] * gpg_13[k];

        t_19[k] = f_4 * pc_y[k] * gpg_14[k];

        t_20[k] = f_2 * gpf0_9[k]
                  - f_3 * gpf1_9[k]
                  + f_4 * pc_z[k] * gpg_14[k];

        t_21[k] = pb_y[k] * gsh0_0[k]
                  - f_9 * pc_y[k] * gsh1_0[k];

        t_22[k] = f_1 * gsg_0[k]
                  + f_4 * pc_y[k] * gpg_15[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pb_y, pc_y, pc_z, gsh0_3, gsh0_5, gsg_1, \
                         gsg_2, gsh1_3, gsh1_5, gpg_15, gpg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_4 * pc_z[k] * gpg_15[k];

        t_24[k] = pb_y[k] * gsh0_3[k]
                  + f_10 * gsg_1[k]
                  - f_9 * pc_y[k] * gsh1_3[k];

        t_25[k] = f_1 * gsg_2[k]
                  + f_4 * pc_y[k] * gpg_17[k];

        t_26[k] = pb_y[k] * gsh0_5[k]
                  - f_9 * pc_y[k] * gsh1_5[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pb_y, pc_y, pc_z, gsh0_6, gsh0_9, gsg_3, \
                         gsg_5, gsh1_6, gsh1_9, gpg_18, gpg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pb_y[k] * gsh0_6[k]
                  + f_11 * gsg_3[k]
                  - f_9 * pc_y[k] * gsh1_6[k];

        t_28[k] = f_4 * pc_z[k] * gpg_18[k];

        t_29[k] = f_1 * gsg_5[k]
                  + f_4 * pc_y[k] * gpg_20[k];

        t_30[k] = pb_y[k] * gsh0_9[k]
                  - f_9 * pc_y[k] * gsh1_9[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pc_x, pc_y, pc_z, fpg_25, fpg_27, gsg_9, \
                         gpg_21, gpg_24, gpg_25, gpg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_0 * fpg_25[k]
                  + f_4 * pc_x[k] * gpg_25[k];

        t_32[k] = f_4 * pc_z[k] * gpg_21[k];

        t_33[k] = f_0 * fpg_27[k]
                  + f_4 * pc_x[k] * gpg_27[k];

        t_34[k] = f_1 * gsg_9[k]
                  + f_4 * pc_y[k] * gpg_24[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pb_y, pc_x, pc_y, pc_z, fpg_29, gsh0_15, gsg_10, \
                         gsh1_15, gpg_25, gpg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_0 * fpg_29[k]
                  + f_4 * pc_x[k] * gpg_29[k];

        t_36[k] = pb_y[k] * gsh0_15[k]
                  + f_12 * gsg_10[k]
                  - f_9 * pc_y[k] * gsh1_15[k];

        t_37[k] = f_4 * pc_z[k] * gpg_25[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_y, pc_y, gsh0_17, gsh0_18, gsh0_20, \
                         gsg_12, gsg_13, gsg_14, gsh1_17, gsh1_18, gsh1_20, \
                         gpg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pb_y[k] * gsh0_17[k]
                  + f_11 * gsg_12[k]
                  - f_9 * pc_y[k] * gsh1_17[k];

        t_39[k] = pb_y[k] * gsh0_18[k]
                  + f_10 * gsg_13[k]
                  - f_9 * pc_y[k] * gsh1_18[k];

        t_40[k] = f_1 * gsg_14[k]
                  + f_4 * pc_y[k] * gpg_29[k];

        t_41[k] = pb_y[k] * gsh0_20[k]
                  - f_9 * pc_y[k] * gsh1_20[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pb_z, pc_y, pc_z, gsh0_0, gsh0_3, \
                         gsg_0, gsh1_0, gsh1_3, gpg_30, gpg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_z[k] * gsh0_0[k]
                  - f_9 * pc_z[k] * gsh1_0[k];

        t_43[k] = f_4 * pc_y[k] * gpg_30[k];

        t_44[k] = f_1 * gsg_0[k]
                  + f_4 * pc_z[k] * gpg_30[k];

        t_45[k] = pb_z[k] * gsh0_3[k]
                  - f_9 * pc_z[k] * gsh1_3[k];

        t_46[k] = f_4 * pc_y[k] * gpg_32[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pb_z, pc_y, pc_z, gsh0_5, gsh0_6, gsg_2, \
                         gsg_3, gsh1_5, gsh1_6, gpg_33, gpg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = pb_z[k] * gsh0_5[k]
                  + f_10 * gsg_2[k]
                  - f_9 * pc_z[k] * gsh1_5[k];

        t_48[k] = pb_z[k] * gsh0_6[k]
                  - f_9 * pc_z[k] * gsh1_6[k];

        t_49[k] = f_1 * gsg_3[k]
                  + f_4 * pc_z[k] * gpg_33[k];

        t_50[k] = f_4 * pc_y[k] * gpg_35[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pb_z, pc_x, pc_z, fpg_40, fpg_42, gsh0_9, \
                         gsg_5, gsg_6, gsh1_9, gpg_36, gpg_40, gpg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pb_z[k] * gsh0_9[k]
                  + f_11 * gsg_5[k]
                  - f_9 * pc_z[k] * gsh1_9[k];

        t_52[k] = f_0 * fpg_40[k]
                  + f_4 * pc_x[k] * gpg_40[k];

        t_53[k] = f_1 * gsg_6[k]
                  + f_4 * pc_z[k] * gpg_36[k];

        t_54[k] = f_0 * fpg_42[k]
                  + f_4 * pc_x[k] * gpg_42[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pb_z, pc_x, pc_y, pc_z, fpg_44, gsh0_15, \
                         gsg_10, gsh1_15, gpg_39, gpg_40, gpg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_4 * pc_y[k] * gpg_39[k];

        t_56[k] = f_0 * fpg_44[k]
                  + f_4 * pc_x[k] * gpg_44[k];

        t_57[k] = pb_z[k] * gsh0_15[k]
                  - f_9 * pc_z[k] * gsh1_15[k];

        t_58[k] = f_1 * gsg_10[k]
                  + f_4 * pc_z[k] * gpg_40[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, pc_y, gpf0_28, gpf0_29, gpf1_28, gpf1_29, gpg_42, \
                         gpg_43, gpg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_7 * gpf0_28[k]
                  - f_8 * gpf1_28[k]
                  + f_4 * pc_y[k] * gpg_42[k];

        t_60[k] = f_5 * gpf0_29[k]
                  - f_6 * gpf1_29[k]
                  + f_4 * pc_y[k] * gpg_43[k];

        t_61[k] = f_4 * pc_y[k] * gpg_44[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_y, pb_z, pc_y, pc_z, fph0_0, fpg_0, \
                         fph1_0, gsh0_20, gsg_14, gsh1_20, gpg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pb_z[k] * gsh0_20[k]
                  + f_12 * gsg_14[k]
                  - f_9 * pc_z[k] * gsh1_20[k];

        t_63[k] = pa_y[k] * fph0_0[k]
                  - f_9 * pc_y[k] * fph1_0[k];

        t_64[k] = f_1 * fpg_0[k]
                  + f_4 * pc_y[k] * gpg_45[k];

        t_65[k] = f_4 * pc_z[k] * gpg_45[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pa_y, pc_y, pc_z, fph0_3, fph0_5, fph0_6, \
                         fpg_1, fpg_3, fph1_3, fph1_5, fph1_6, gpg_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_y[k] * fph0_3[k]
                  + f_10 * fpg_1[k]
                  - f_9 * pc_y[k] * fph1_3[k];

        t_67[k] = f_4 * pc_z[k] * gpg_46[k];

        t_68[k] = pa_y[k] * fph0_5[k]
                  - f_9 * pc_y[k] * fph1_5[k];

        t_69[k] = pa_y[k] * fph0_6[k]
                  + f_11 * fpg_3[k]
                  - f_9 * pc_y[k] * fph1_6[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pa_y, pc_x, pc_y, pc_z, fph0_9, fpg_5, \
                         fpg_55, fph1_9, gsg_25, gpg_48, gpg_50, \
                         gpg_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_4 * pc_z[k] * gpg_48[k];

        t_71[k] = f_1 * fpg_5[k]
                  + f_4 * pc_y[k] * gpg_50[k];

        t_72[k] = pa_y[k] * fph0_9[k]
                  - f_9 * pc_y[k] * fph1_9[k];

        t_73[k] = f_11 * fpg_55[k]
                  + f_1 * gsg_25[k]
                  + f_4 * pc_x[k] * gpg_55[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, pc_x, pc_z, fpg_57, fpg_58, gsg_27, gsg_28, gpg_51, \
                         gpg_57, gpg_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_4 * pc_z[k] * gpg_51[k];

        t_75[k] = f_11 * fpg_57[k]
                  + f_1 * gsg_27[k]
                  + f_4 * pc_x[k] * gpg_57[k];

        t_76[k] = f_11 * fpg_58[k]
                  + f_1 * gsg_28[k]
                  + f_4 * pc_x[k] * gpg_58[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_y, pc_y, pc_z, fph0_14, fpg_10, fph1_14, \
                         gpf0_36, gpf1_36, gpg_55, gpg_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pa_y[k] * fph0_14[k]
                  - f_9 * pc_y[k] * fph1_14[k];

        t_78[k] = f_1 * fpg_10[k]
                  + f_2 * gpf0_36[k]
                  - f_3 * gpf1_36[k]
                  + f_4 * pc_y[k] * gpg_55[k];

        t_79[k] = f_4 * pc_z[k] * gpg_55[k];

        t_80[k] = f_5 * gpf0_36[k]
                  - f_6 * gpf1_36[k]
                  + f_4 * pc_z[k] * gpg_56[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pa_y, pc_y, pc_z, fph0_20, fpg_14, fph1_20, \
                         gpf0_37, gpf1_37, gpg_57, gpg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_7 * gpf0_37[k]
                  - f_8 * gpf1_37[k]
                  + f_4 * pc_z[k] * gpg_57[k];

        t_82[k] = f_1 * fpg_14[k]
                  + f_4 * pc_y[k] * gpg_59[k];

        t_83[k] = pa_y[k] * fph0_20[k]
                  - f_9 * pc_y[k] * fph1_20[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pc_x, pc_y, pc_z, fpg_15, fpg_60, gsg_15, gpf0_40, \
                         gpf1_40, gpg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_11 * fpg_60[k]
                  + f_2 * gpf0_40[k]
                  - f_3 * gpf1_40[k]
                  + f_4 * pc_x[k] * gpg_60[k];

        t_85[k] = f_1 * fpg_15[k]
                  + f_1 * gsg_15[k]
                  + f_4 * pc_y[k] * gpg_60[k];

        t_86[k] = f_4 * pc_z[k] * gpg_60[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pc_x, pc_z, fpg_63, gpf0_40, gpf0_43, gpf1_40, \
                         gpf1_43, gpg_61, gpg_62, gpg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_11 * fpg_63[k]
                  + f_7 * gpf0_43[k]
                  - f_8 * gpf1_43[k]
                  + f_4 * pc_x[k] * gpg_63[k];

        t_88[k] = f_4 * pc_z[k] * gpg_61[k];

        t_89[k] = f_5 * gpf0_40[k]
                  - f_6 * gpf1_40[k]
                  + f_4 * pc_z[k] * gpg_62[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, pc_x, pc_y, pc_z, fpg_20, fpg_66, gsg_20, gpf0_46, \
                         gpf1_46, gpg_63, gpg_65, gpg_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_11 * fpg_66[k]
                  + f_5 * gpf0_46[k]
                  - f_6 * gpf1_46[k]
                  + f_4 * pc_x[k] * gpg_66[k];

        t_91[k] = f_4 * pc_z[k] * gpg_63[k];

        t_92[k] = f_1 * fpg_20[k]
                  + f_1 * gsg_20[k]
                  + f_4 * pc_y[k] * gpg_65[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pc_x, pc_z, fpg_70, fpg_72, gpf0_42, gpf1_42, \
                         gpg_65, gpg_66, gpg_70, gpg_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_7 * gpf0_42[k]
                  - f_8 * gpf1_42[k]
                  + f_4 * pc_z[k] * gpg_65[k];

        t_94[k] = f_11 * fpg_70[k]
                  + f_4 * pc_x[k] * gpg_70[k];

        t_95[k] = f_4 * pc_z[k] * gpg_66[k];

        t_96[k] = f_11 * fpg_72[k]
                  + f_4 * pc_x[k] * gpg_72[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pa_x, pc_x, pc_z, dph0_99, dph1_99, fph0_99, \
                         fpg_73, fpg_74, fph1_99, gpg_70, gpg_73, \
                         gpg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_11 * fpg_73[k]
                  + f_4 * pc_x[k] * gpg_73[k];

        t_98[k] = f_11 * fpg_74[k]
                  + f_4 * pc_x[k] * gpg_74[k];

        t_99[k] = f_13 * dph0_99[k]
                  - f_14 * dph1_99[k]
                  + pa_x[k] * fph0_99[k]
                  - f_9 * pc_x[k] * fph1_99[k];

        t_100[k] = f_4 * pc_z[k] * gpg_70[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, pc_y, pc_z, fpg_29, gsg_29, gpf0_46, gpf0_47, \
                         gpf1_46, gpf1_47, gpg_71, gpg_72, gpg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_5 * gpf0_46[k]
                   - f_6 * gpf1_46[k]
                   + f_4 * pc_z[k] * gpg_71[k];

        t_102[k] = f_7 * gpf0_47[k]
                   - f_8 * gpf1_47[k]
                   + f_4 * pc_z[k] * gpg_72[k];

        t_103[k] = f_1 * fpg_29[k]
                   + f_1 * gsg_29[k]
                   + f_4 * pc_y[k] * gpg_74[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pa_y, pc_y, pc_z, fph0_42, fpg_30, \
                         fph1_42, gsg_15, gpf0_49, gpf1_49, gpg_74, \
                         gpg_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_2 * gpf0_49[k]
                   - f_3 * gpf1_49[k]
                   + f_4 * pc_z[k] * gpg_74[k];

        t_105[k] = pa_y[k] * fph0_42[k]
                   - f_9 * pc_y[k] * fph1_42[k];

        t_106[k] = f_1 * fpg_30[k]
                   + f_4 * pc_y[k] * gpg_75[k];

        t_107[k] = f_1 * gsg_15[k]
                   + f_4 * pc_z[k] * gpg_75[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, pa_y, pb_z, pc_y, pc_z, fph0_47, fph1_47, \
                         gsh0_24, gsh0_27, gsg_16, gsh1_24, gsh1_27, \
                         gpg_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = pb_z[k] * gsh0_24[k]
                   - f_9 * pc_z[k] * gsh1_24[k];

        t_109[k] = f_1 * gsg_16[k]
                   + f_4 * pc_z[k] * gpg_76[k];

        t_110[k] = pa_y[k] * fph0_47[k]
                   - f_9 * pc_y[k] * fph1_47[k];

        t_111[k] = pb_z[k] * gsh0_27[k]
                   - f_9 * pc_z[k] * gsh1_27[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, pa_y, pc_x, pc_y, pc_z, fph0_51, fpg_35, \
                         fpg_85, fph1_51, gsg_18, gpg_78, gpg_80, \
                         gpg_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = f_1 * gsg_18[k]
                   + f_4 * pc_z[k] * gpg_78[k];

        t_113[k] = f_1 * fpg_35[k]
                   + f_4 * pc_y[k] * gpg_80[k];

        t_114[k] = pa_y[k] * fph0_51[k]
                   - f_9 * pc_y[k] * fph1_51[k];

        t_115[k] = f_11 * fpg_85[k]
                   + f_4 * pc_x[k] * gpg_85[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pc_x, pc_z, fpg_87, fpg_88, fpg_89, \
                         gsg_21, gpg_81, gpg_87, gpg_88, gpg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = f_1 * gsg_21[k]
                   + f_4 * pc_z[k] * gpg_81[k];

        t_117[k] = f_11 * fpg_87[k]
                   + f_4 * pc_x[k] * gpg_87[k];

        t_118[k] = f_11 * fpg_88[k]
                   + f_4 * pc_x[k] * gpg_88[k];

        t_119[k] = f_11 * fpg_89[k]
                   + f_4 * pc_x[k] * gpg_89[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pb_z, pc_z, gsh0_36, gsh0_38, gsh0_39, \
                         gsg_25, gsg_26, gsg_27, gsh1_36, gsh1_38, gsh1_39, \
                         gpg_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = pb_z[k] * gsh0_36[k]
                   - f_9 * pc_z[k] * gsh1_36[k];

        t_121[k] = f_1 * gsg_25[k]
                   + f_4 * pc_z[k] * gpg_85[k];

        t_122[k] = pb_z[k] * gsh0_38[k]
                   + f_10 * gsg_26[k]
                   - f_9 * pc_z[k] * gsh1_38[k];

        t_123[k] = pb_z[k] * gsh0_39[k]
                   + f_11 * gsg_27[k]
                   - f_9 * pc_z[k] * gsh1_39[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pa_y, pa_z, pc_y, pc_z, fph0_0, fph0_62, \
                         fpg_44, fph1_0, fph1_62, gpg_89, gpg_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_1 * fpg_44[k]
                   + f_4 * pc_y[k] * gpg_89[k];

        t_125[k] = pa_y[k] * fph0_62[k]
                   - f_9 * pc_y[k] * fph1_62[k];

        t_126[k] = pa_z[k] * fph0_0[k]
                   - f_9 * pc_z[k] * fph1_0[k];

        t_127[k] = f_4 * pc_y[k] * gpg_90[k];
    }
}

static auto
compute_prim_gph_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t dph0, const size_t dph1,
                                                          const size_t fph0, const size_t fpg,
                                                          const size_t fph1, const size_t gsh0,
                                                          const size_t gsg, const size_t gsh1,
                                                          const size_t gpf0, const size_t gpf1,
                                                          const size_t gpg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 2.0 / gamma;
    const auto f_3 = 2.0 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = gamma / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_13 = 1.0 / p;
    const auto f_14 = gamma / (p * q);
    const auto f_15 = 1.5 / gamma;
    const auto f_16 = 1.5 * p / (gamma * q);
    const auto f_17 = 0.5 / p;
    const auto f_18 = 0.5 * gamma / (p * q);

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
    auto *t_246 = buffer.data(target + 246);
    auto *t_247 = buffer.data(target + 247);
    auto *t_248 = buffer.data(target + 248);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dph0_0 = buffer.data(dph0 + 0);
    const auto *dph0_188 = buffer.data(dph0 + 188);
    const auto *dph0_225 = buffer.data(dph0 + 225);

    const auto *dph1_0 = buffer.data(dph1 + 0);
    const auto *dph1_188 = buffer.data(dph1 + 188);
    const auto *dph1_225 = buffer.data(dph1 + 225);

    const auto *fph0_3 = buffer.data(fph0 + 3);
    const auto *fph0_5 = buffer.data(fph0 + 5);
    const auto *fph0_6 = buffer.data(fph0 + 6);
    const auto *fph0_9 = buffer.data(fph0 + 9);
    const auto *fph0_10 = buffer.data(fph0 + 10);
    const auto *fph0_15 = buffer.data(fph0 + 15);
    const auto *fph0_21 = buffer.data(fph0 + 21);
    const auto *fph0_24 = buffer.data(fph0 + 24);
    const auto *fph0_27 = buffer.data(fph0 + 27);
    const auto *fph0_28 = buffer.data(fph0 + 28);
    const auto *fph0_36 = buffer.data(fph0 + 36);
    const auto *fph0_63 = buffer.data(fph0 + 63);
    const auto *fph0_188 = buffer.data(fph0 + 188);
    const auto *fph0_225 = buffer.data(fph0 + 225);

    const auto *fpg_0 = buffer.data(fpg + 0);
    const auto *fpg_2 = buffer.data(fpg + 2);
    const auto *fpg_5 = buffer.data(fpg + 5);
    const auto *fpg_14 = buffer.data(fpg + 14);
    const auto *fpg_15 = buffer.data(fpg + 15);
    const auto *fpg_18 = buffer.data(fpg + 18);
    const auto *fpg_30 = buffer.data(fpg + 30);
    const auto *fpg_45 = buffer.data(fpg + 45);
    const auto *fpg_50 = buffer.data(fpg + 50);
    const auto *fpg_55 = buffer.data(fpg + 55);
    const auto *fpg_59 = buffer.data(fpg + 59);
    const auto *fpg_60 = buffer.data(fpg + 60);
    const auto *fpg_65 = buffer.data(fpg + 65);
    const auto *fpg_74 = buffer.data(fpg + 74);
    const auto *fpg_75 = buffer.data(fpg + 75);
    const auto *fpg_80 = buffer.data(fpg + 80);
    const auto *fpg_101 = buffer.data(fpg + 101);
    const auto *fpg_102 = buffer.data(fpg + 102);
    const auto *fpg_104 = buffer.data(fpg + 104);
    const auto *fpg_115 = buffer.data(fpg + 115);
    const auto *fpg_116 = buffer.data(fpg + 116);
    const auto *fpg_117 = buffer.data(fpg + 117);
    const auto *fpg_119 = buffer.data(fpg + 119);
    const auto *fpg_120 = buffer.data(fpg + 120);
    const auto *fpg_125 = buffer.data(fpg + 125);
    const auto *fpg_129 = buffer.data(fpg + 129);
    const auto *fpg_130 = buffer.data(fpg + 130);
    const auto *fpg_131 = buffer.data(fpg + 131);
    const auto *fpg_132 = buffer.data(fpg + 132);
    const auto *fpg_134 = buffer.data(fpg + 134);
    const auto *fpg_138 = buffer.data(fpg + 138);
    const auto *fpg_141 = buffer.data(fpg + 141);
    const auto *fpg_145 = buffer.data(fpg + 145);
    const auto *fpg_147 = buffer.data(fpg + 147);
    const auto *fpg_148 = buffer.data(fpg + 148);
    const auto *fpg_149 = buffer.data(fpg + 149);
    const auto *fpg_150 = buffer.data(fpg + 150);
    const auto *fpg_153 = buffer.data(fpg + 153);
    const auto *fpg_156 = buffer.data(fpg + 156);
    const auto *fpg_160 = buffer.data(fpg + 160);
    const auto *fpg_162 = buffer.data(fpg + 162);
    const auto *fpg_163 = buffer.data(fpg + 163);
    const auto *fpg_164 = buffer.data(fpg + 164);
    const auto *fpg_175 = buffer.data(fpg + 175);
    const auto *fpg_177 = buffer.data(fpg + 177);
    const auto *fpg_178 = buffer.data(fpg + 178);
    const auto *fpg_179 = buffer.data(fpg + 179);

    const auto *fph1_3 = buffer.data(fph1 + 3);
    const auto *fph1_5 = buffer.data(fph1 + 5);
    const auto *fph1_6 = buffer.data(fph1 + 6);
    const auto *fph1_9 = buffer.data(fph1 + 9);
    const auto *fph1_10 = buffer.data(fph1 + 10);
    const auto *fph1_15 = buffer.data(fph1 + 15);
    const auto *fph1_21 = buffer.data(fph1 + 21);
    const auto *fph1_24 = buffer.data(fph1 + 24);
    const auto *fph1_27 = buffer.data(fph1 + 27);
    const auto *fph1_28 = buffer.data(fph1 + 28);
    const auto *fph1_36 = buffer.data(fph1 + 36);
    const auto *fph1_63 = buffer.data(fph1 + 63);
    const auto *fph1_188 = buffer.data(fph1 + 188);
    const auto *fph1_225 = buffer.data(fph1 + 225);

    const auto *gsh0_47 = buffer.data(gsh0 + 47);
    const auto *gsh0_51 = buffer.data(gsh0 + 51);
    const auto *gsh0_58 = buffer.data(gsh0 + 58);
    const auto *gsh0_59 = buffer.data(gsh0 + 59);
    const auto *gsh0_60 = buffer.data(gsh0 + 60);
    const auto *gsh0_62 = buffer.data(gsh0 + 62);
    const auto *gsh0_63 = buffer.data(gsh0 + 63);
    const auto *gsh0_66 = buffer.data(gsh0 + 66);
    const auto *gsh0_68 = buffer.data(gsh0 + 68);
    const auto *gsh0_69 = buffer.data(gsh0 + 69);
    const auto *gsh0_72 = buffer.data(gsh0 + 72);
    const auto *gsh0_78 = buffer.data(gsh0 + 78);
    const auto *gsh0_80 = buffer.data(gsh0 + 80);

    const auto *gsg_30 = buffer.data(gsg + 30);
    const auto *gsg_32 = buffer.data(gsg + 32);
    const auto *gsg_35 = buffer.data(gsg + 35);
    const auto *gsg_39 = buffer.data(gsg + 39);
    const auto *gsg_41 = buffer.data(gsg + 41);
    const auto *gsg_42 = buffer.data(gsg + 42);
    const auto *gsg_43 = buffer.data(gsg + 43);
    const auto *gsg_44 = buffer.data(gsg + 44);
    const auto *gsg_45 = buffer.data(gsg + 45);
    const auto *gsg_46 = buffer.data(gsg + 46);
    const auto *gsg_47 = buffer.data(gsg + 47);
    const auto *gsg_48 = buffer.data(gsg + 48);
    const auto *gsg_50 = buffer.data(gsg + 50);
    const auto *gsg_51 = buffer.data(gsg + 51);
    const auto *gsg_55 = buffer.data(gsg + 55);
    const auto *gsg_56 = buffer.data(gsg + 56);
    const auto *gsg_57 = buffer.data(gsg + 57);
    const auto *gsg_58 = buffer.data(gsg + 58);
    const auto *gsg_59 = buffer.data(gsg + 59);

    const auto *gsh1_47 = buffer.data(gsh1 + 47);
    const auto *gsh1_51 = buffer.data(gsh1 + 51);
    const auto *gsh1_58 = buffer.data(gsh1 + 58);
    const auto *gsh1_59 = buffer.data(gsh1 + 59);
    const auto *gsh1_60 = buffer.data(gsh1 + 60);
    const auto *gsh1_62 = buffer.data(gsh1 + 62);
    const auto *gsh1_63 = buffer.data(gsh1 + 63);
    const auto *gsh1_66 = buffer.data(gsh1 + 66);
    const auto *gsh1_68 = buffer.data(gsh1 + 68);
    const auto *gsh1_69 = buffer.data(gsh1 + 69);
    const auto *gsh1_72 = buffer.data(gsh1 + 72);
    const auto *gsh1_78 = buffer.data(gsh1 + 78);
    const auto *gsh1_80 = buffer.data(gsh1 + 80);

    const auto *gpf0_62 = buffer.data(gpf0 + 62);
    const auto *gpf0_67 = buffer.data(gpf0 + 67);
    const auto *gpf0_68 = buffer.data(gpf0 + 68);
    const auto *gpf0_69 = buffer.data(gpf0 + 69);
    const auto *gpf0_80 = buffer.data(gpf0 + 80);
    const auto *gpf0_81 = buffer.data(gpf0 + 81);
    const auto *gpf0_82 = buffer.data(gpf0 + 82);
    const auto *gpf0_85 = buffer.data(gpf0 + 85);
    const auto *gpf0_86 = buffer.data(gpf0 + 86);
    const auto *gpf0_87 = buffer.data(gpf0 + 87);
    const auto *gpf0_88 = buffer.data(gpf0 + 88);
    const auto *gpf0_89 = buffer.data(gpf0 + 89);
    const auto *gpf0_90 = buffer.data(gpf0 + 90);
    const auto *gpf0_92 = buffer.data(gpf0 + 92);
    const auto *gpf0_93 = buffer.data(gpf0 + 93);
    const auto *gpf0_96 = buffer.data(gpf0 + 96);
    const auto *gpf0_97 = buffer.data(gpf0 + 97);
    const auto *gpf0_99 = buffer.data(gpf0 + 99);
    const auto *gpf0_100 = buffer.data(gpf0 + 100);
    const auto *gpf0_102 = buffer.data(gpf0 + 102);
    const auto *gpf0_103 = buffer.data(gpf0 + 103);
    const auto *gpf0_106 = buffer.data(gpf0 + 106);
    const auto *gpf0_107 = buffer.data(gpf0 + 107);
    const auto *gpf0_109 = buffer.data(gpf0 + 109);

    const auto *gpf1_62 = buffer.data(gpf1 + 62);
    const auto *gpf1_67 = buffer.data(gpf1 + 67);
    const auto *gpf1_68 = buffer.data(gpf1 + 68);
    const auto *gpf1_69 = buffer.data(gpf1 + 69);
    const auto *gpf1_80 = buffer.data(gpf1 + 80);
    const auto *gpf1_81 = buffer.data(gpf1 + 81);
    const auto *gpf1_82 = buffer.data(gpf1 + 82);
    const auto *gpf1_85 = buffer.data(gpf1 + 85);
    const auto *gpf1_86 = buffer.data(gpf1 + 86);
    const auto *gpf1_87 = buffer.data(gpf1 + 87);
    const auto *gpf1_88 = buffer.data(gpf1 + 88);
    const auto *gpf1_89 = buffer.data(gpf1 + 89);
    const auto *gpf1_90 = buffer.data(gpf1 + 90);
    const auto *gpf1_92 = buffer.data(gpf1 + 92);
    const auto *gpf1_93 = buffer.data(gpf1 + 93);
    const auto *gpf1_96 = buffer.data(gpf1 + 96);
    const auto *gpf1_97 = buffer.data(gpf1 + 97);
    const auto *gpf1_99 = buffer.data(gpf1 + 99);
    const auto *gpf1_100 = buffer.data(gpf1 + 100);
    const auto *gpf1_102 = buffer.data(gpf1 + 102);
    const auto *gpf1_103 = buffer.data(gpf1 + 103);
    const auto *gpf1_106 = buffer.data(gpf1 + 106);
    const auto *gpf1_107 = buffer.data(gpf1 + 107);
    const auto *gpf1_109 = buffer.data(gpf1 + 109);

    const auto *gpg_90 = buffer.data(gpg + 90);
    const auto *gpg_92 = buffer.data(gpg + 92);
    const auto *gpg_94 = buffer.data(gpg + 94);
    const auto *gpg_95 = buffer.data(gpg + 95);
    const auto *gpg_99 = buffer.data(gpg + 99);
    const auto *gpg_101 = buffer.data(gpg + 101);
    const auto *gpg_102 = buffer.data(gpg + 102);
    const auto *gpg_103 = buffer.data(gpg + 103);
    const auto *gpg_104 = buffer.data(gpg + 104);
    const auto *gpg_105 = buffer.data(gpg + 105);
    const auto *gpg_107 = buffer.data(gpg + 107);
    const auto *gpg_110 = buffer.data(gpg + 110);
    const auto *gpg_114 = buffer.data(gpg + 114);
    const auto *gpg_115 = buffer.data(gpg + 115);
    const auto *gpg_116 = buffer.data(gpg + 116);
    const auto *gpg_117 = buffer.data(gpg + 117);
    const auto *gpg_119 = buffer.data(gpg + 119);
    const auto *gpg_120 = buffer.data(gpg + 120);
    const auto *gpg_121 = buffer.data(gpg + 121);
    const auto *gpg_122 = buffer.data(gpg + 122);
    const auto *gpg_123 = buffer.data(gpg + 123);
    const auto *gpg_124 = buffer.data(gpg + 124);
    const auto *gpg_125 = buffer.data(gpg + 125);
    const auto *gpg_129 = buffer.data(gpg + 129);
    const auto *gpg_130 = buffer.data(gpg + 130);
    const auto *gpg_131 = buffer.data(gpg + 131);
    const auto *gpg_132 = buffer.data(gpg + 132);
    const auto *gpg_133 = buffer.data(gpg + 133);
    const auto *gpg_134 = buffer.data(gpg + 134);
    const auto *gpg_135 = buffer.data(gpg + 135);
    const auto *gpg_136 = buffer.data(gpg + 136);
    const auto *gpg_137 = buffer.data(gpg + 137);
    const auto *gpg_138 = buffer.data(gpg + 138);
    const auto *gpg_140 = buffer.data(gpg + 140);
    const auto *gpg_141 = buffer.data(gpg + 141);
    const auto *gpg_145 = buffer.data(gpg + 145);
    const auto *gpg_146 = buffer.data(gpg + 146);
    const auto *gpg_147 = buffer.data(gpg + 147);
    const auto *gpg_148 = buffer.data(gpg + 148);
    const auto *gpg_149 = buffer.data(gpg + 149);
    const auto *gpg_150 = buffer.data(gpg + 150);
    const auto *gpg_151 = buffer.data(gpg + 151);
    const auto *gpg_152 = buffer.data(gpg + 152);
    const auto *gpg_153 = buffer.data(gpg + 153);
    const auto *gpg_155 = buffer.data(gpg + 155);
    const auto *gpg_156 = buffer.data(gpg + 156);
    const auto *gpg_160 = buffer.data(gpg + 160);
    const auto *gpg_161 = buffer.data(gpg + 161);
    const auto *gpg_162 = buffer.data(gpg + 162);
    const auto *gpg_163 = buffer.data(gpg + 163);
    const auto *gpg_164 = buffer.data(gpg + 164);
    const auto *gpg_165 = buffer.data(gpg + 165);
    const auto *gpg_166 = buffer.data(gpg + 166);
    const auto *gpg_168 = buffer.data(gpg + 168);
    const auto *gpg_170 = buffer.data(gpg + 170);
    const auto *gpg_171 = buffer.data(gpg + 171);
    const auto *gpg_175 = buffer.data(gpg + 175);
    const auto *gpg_177 = buffer.data(gpg + 177);
    const auto *gpg_178 = buffer.data(gpg + 178);
    const auto *gpg_179 = buffer.data(gpg + 179);

#pragma omp simd aligned(t_128, t_129, t_130, t_131, pa_z, pc_y, pc_z, fph0_3, fph0_5, fpg_0, \
                         fpg_2, fph1_3, fph1_5, gpg_90, gpg_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_1 * fpg_0[k]
                   + f_4 * pc_z[k] * gpg_90[k];

        t_129[k] = pa_z[k] * fph0_3[k]
                   - f_9 * pc_z[k] * fph1_3[k];

        t_130[k] = f_4 * pc_y[k] * gpg_92[k];

        t_131[k] = pa_z[k] * fph0_5[k]
                   + f_10 * fpg_2[k]
                   - f_9 * pc_z[k] * fph1_5[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pa_z, pc_y, pc_z, fph0_6, fph0_9, fpg_5, \
                         fph1_6, fph1_9, gpf0_62, gpf1_62, gpg_94, \
                         gpg_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = pa_z[k] * fph0_6[k]
                   - f_9 * pc_z[k] * fph1_6[k];

        t_133[k] = f_5 * gpf0_62[k]
                   - f_6 * gpf1_62[k]
                   + f_4 * pc_y[k] * gpg_94[k];

        t_134[k] = f_4 * pc_y[k] * gpg_95[k];

        t_135[k] = pa_z[k] * fph0_9[k]
                   + f_11 * fpg_5[k]
                   - f_9 * pc_z[k] * fph1_9[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pa_z, pc_x, pc_z, fph0_10, fpg_101, fpg_102, \
                         fph1_10, gsg_41, gsg_42, gpg_101, gpg_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = pa_z[k] * fph0_10[k]
                   - f_9 * pc_z[k] * fph1_10[k];

        t_137[k] = f_11 * fpg_101[k]
                   + f_1 * gsg_41[k]
                   + f_4 * pc_x[k] * gpg_101[k];

        t_138[k] = f_11 * fpg_102[k]
                   + f_1 * gsg_42[k]
                   + f_4 * pc_x[k] * gpg_102[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, pa_z, pc_x, pc_y, pc_z, fph0_15, fpg_104, \
                         fph1_15, gsg_44, gpg_99, gpg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_4 * pc_y[k] * gpg_99[k];

        t_140[k] = f_11 * fpg_104[k]
                   + f_1 * gsg_44[k]
                   + f_4 * pc_x[k] * gpg_104[k];

        t_141[k] = pa_z[k] * fph0_15[k]
                   - f_9 * pc_z[k] * fph1_15[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pc_y, gpf0_67, gpf0_68, gpf0_69, gpf1_67, \
                         gpf1_68, gpf1_69, gpg_101, gpg_102, gpg_103, \
                         gpg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_15 * gpf0_67[k]
                   - f_16 * gpf1_67[k]
                   + f_4 * pc_y[k] * gpg_101[k];

        t_143[k] = f_7 * gpf0_68[k]
                   - f_8 * gpf1_68[k]
                   + f_4 * pc_y[k] * gpg_102[k];

        t_144[k] = f_5 * gpf0_69[k]
                   - f_6 * gpf1_69[k]
                   + f_4 * pc_y[k] * gpg_103[k];

        t_145[k] = f_4 * pc_y[k] * gpg_104[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pa_z, pc_y, pc_z, fph0_21, fpg_14, \
                         fpg_15, fph1_21, gsg_30, gpf0_69, gpf1_69, gpg_104, \
                         gpg_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_1 * fpg_14[k]
                   + f_2 * gpf0_69[k]
                   - f_3 * gpf1_69[k]
                   + f_4 * pc_z[k] * gpg_104[k];

        t_147[k] = pa_z[k] * fph0_21[k]
                   - f_9 * pc_z[k] * fph1_21[k];

        t_148[k] = f_1 * gsg_30[k]
                   + f_4 * pc_y[k] * gpg_105[k];

        t_149[k] = f_1 * fpg_15[k]
                   + f_4 * pc_z[k] * gpg_105[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, t_153, pa_z, pb_y, pc_y, pc_z, fph0_24, fph0_27, \
                         fph1_24, fph1_27, gsh0_47, gsg_32, gsh1_47, \
                         gpg_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pa_z[k] * fph0_24[k]
                   - f_9 * pc_z[k] * fph1_24[k];

        t_151[k] = f_1 * gsg_32[k]
                   + f_4 * pc_y[k] * gpg_107[k];

        t_152[k] = pb_y[k] * gsh0_47[k]
                   - f_9 * pc_y[k] * gsh1_47[k];

        t_153[k] = pa_z[k] * fph0_27[k]
                   - f_9 * pc_z[k] * fph1_27[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, pa_z, pb_y, pc_y, pc_z, fph0_28, fpg_18, \
                         fph1_28, gsh0_51, gsg_35, gsh1_51, gpg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = pa_z[k] * fph0_28[k]
                   + f_1 * fpg_18[k]
                   - f_9 * pc_z[k] * fph1_28[k];

        t_155[k] = f_1 * gsg_35[k]
                   + f_4 * pc_y[k] * gpg_110[k];

        t_156[k] = pb_y[k] * gsh0_51[k]
                   - f_9 * pc_y[k] * gsh1_51[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, pc_x, pc_y, fpg_115, fpg_116, fpg_117, \
                         gsg_39, gpg_114, gpg_115, gpg_116, gpg_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_11 * fpg_115[k]
                   + f_4 * pc_x[k] * gpg_115[k];

        t_158[k] = f_11 * fpg_116[k]
                   + f_4 * pc_x[k] * gpg_116[k];

        t_159[k] = f_11 * fpg_117[k]
                   + f_4 * pc_x[k] * gpg_117[k];

        t_160[k] = f_1 * gsg_39[k]
                   + f_4 * pc_y[k] * gpg_114[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, pa_z, pb_y, pc_x, pc_y, pc_z, fph0_36, fpg_119, \
                         fph1_36, gsh0_58, gsg_41, gsh1_58, gpg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_11 * fpg_119[k]
                   + f_4 * pc_x[k] * gpg_119[k];

        t_162[k] = pa_z[k] * fph0_36[k]
                   - f_9 * pc_z[k] * fph1_36[k];

        t_163[k] = pb_y[k] * gsh0_58[k]
                   + f_0 * gsg_41[k]
                   - f_9 * pc_y[k] * gsh1_58[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pb_y, pc_y, gsh0_59, gsh0_60, gsh0_62, \
                         gsg_42, gsg_43, gsg_44, gsh1_59, gsh1_60, gsh1_62, \
                         gpg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = pb_y[k] * gsh0_59[k]
                   + f_11 * gsg_42[k]
                   - f_9 * pc_y[k] * gsh1_59[k];

        t_165[k] = pb_y[k] * gsh0_60[k]
                   + f_10 * gsg_43[k]
                   - f_9 * pc_y[k] * gsh1_60[k];

        t_166[k] = f_1 * gsg_44[k]
                   + f_4 * pc_y[k] * gpg_119[k];

        t_167[k] = pb_y[k] * gsh0_62[k]
                   - f_9 * pc_y[k] * gsh1_62[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, t_172, pc_x, pc_y, pc_z, fpg_30, fpg_120, \
                         gsg_30, gpf0_80, gpf1_80, gpg_120, gpg_121, \
                         gpg_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_11 * fpg_120[k]
                   + f_2 * gpf0_80[k]
                   - f_3 * gpf1_80[k]
                   + f_4 * pc_x[k] * gpg_120[k];

        t_169[k] = f_4 * pc_y[k] * gpg_120[k];

        t_170[k] = f_1 * fpg_30[k]
                   + f_1 * gsg_30[k]
                   + f_4 * pc_z[k] * gpg_120[k];

        t_171[k] = f_5 * gpf0_80[k]
                   - f_6 * gpf1_80[k]
                   + f_4 * pc_y[k] * gpg_121[k];

        t_172[k] = f_4 * pc_y[k] * gpg_122[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, pc_x, pc_y, fpg_125, gpf0_81, gpf0_82, \
                         gpf0_85, gpf1_81, gpf1_82, gpf1_85, gpg_123, gpg_124, \
                         gpg_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_11 * fpg_125[k]
                   + f_7 * gpf0_85[k]
                   - f_8 * gpf1_85[k]
                   + f_4 * pc_x[k] * gpg_125[k];

        t_174[k] = f_7 * gpf0_81[k]
                   - f_8 * gpf1_81[k]
                   + f_4 * pc_y[k] * gpg_123[k];

        t_175[k] = f_5 * gpf0_82[k]
                   - f_6 * gpf1_82[k]
                   + f_4 * pc_y[k] * gpg_124[k];

        t_176[k] = f_4 * pc_y[k] * gpg_125[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, pc_x, fpg_129, fpg_130, fpg_131, fpg_132, \
                         gpf0_89, gpf1_89, gpg_129, gpg_130, gpg_131, \
                         gpg_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_11 * fpg_129[k]
                   + f_5 * gpf0_89[k]
                   - f_6 * gpf1_89[k]
                   + f_4 * pc_x[k] * gpg_129[k];

        t_178[k] = f_11 * fpg_130[k]
                   + f_4 * pc_x[k] * gpg_130[k];

        t_179[k] = f_11 * fpg_131[k]
                   + f_4 * pc_x[k] * gpg_131[k];

        t_180[k] = f_11 * fpg_132[k]
                   + f_4 * pc_x[k] * gpg_132[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pc_x, pc_y, fpg_134, gpf0_86, gpf0_87, \
                         gpf1_86, gpf1_87, gpg_129, gpg_130, gpg_131, \
                         gpg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_4 * pc_y[k] * gpg_129[k];

        t_182[k] = f_11 * fpg_134[k]
                   + f_4 * pc_x[k] * gpg_134[k];

        t_183[k] = f_2 * gpf0_86[k]
                   - f_3 * gpf1_86[k]
                   + f_4 * pc_y[k] * gpg_130[k];

        t_184[k] = f_15 * gpf0_87[k]
                   - f_16 * gpf1_87[k]
                   + f_4 * pc_y[k] * gpg_131[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, pc_y, gpf0_88, gpf0_89, gpf1_88, gpf1_89, \
                         gpg_132, gpg_133, gpg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_7 * gpf0_88[k]
                   - f_8 * gpf1_88[k]
                   + f_4 * pc_y[k] * gpg_132[k];

        t_186[k] = f_5 * gpf0_89[k]
                   - f_6 * gpf1_89[k]
                   + f_4 * pc_y[k] * gpg_133[k];

        t_187[k] = f_4 * pc_y[k] * gpg_134[k];
    }

#pragma omp simd aligned(t_188, t_189, pa_x, pa_y, pc_x, pc_y, dph0_0, dph0_188, dph1_0, \
                         dph1_188, fph0_63, fph0_188, fph1_63, \
                         fph1_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_13 * dph0_188[k]
                   - f_14 * dph1_188[k]
                   + pa_x[k] * fph0_188[k]
                   - f_9 * pc_x[k] * fph1_188[k];

        t_189[k] = f_17 * dph0_0[k]
                   - f_18 * dph1_0[k]
                   + pa_y[k] * fph0_63[k]
                   - f_9 * pc_y[k] * fph1_63[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, pc_x, pc_y, pc_z, fpg_45, fpg_138, \
                         gsg_48, gpf0_93, gpf1_93, gpg_135, gpg_136, \
                         gpg_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = f_10 * fpg_45[k]
                   + f_4 * pc_y[k] * gpg_135[k];

        t_191[k] = f_4 * pc_z[k] * gpg_135[k];

        t_192[k] = f_10 * fpg_138[k]
                   + f_1 * gsg_48[k]
                   + f_7 * gpf0_93[k]
                   - f_8 * gpf1_93[k]
                   + f_4 * pc_x[k] * gpg_138[k];

        t_193[k] = f_4 * pc_z[k] * gpg_136[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, pc_x, pc_z, fpg_141, gsg_51, gpf0_90, gpf0_96, \
                         gpf1_90, gpf1_96, gpg_137, gpg_138, gpg_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = f_5 * gpf0_90[k]
                   - f_6 * gpf1_90[k]
                   + f_4 * pc_z[k] * gpg_137[k];

        t_195[k] = f_10 * fpg_141[k]
                   + f_1 * gsg_51[k]
                   + f_5 * gpf0_96[k]
                   - f_6 * gpf1_96[k]
                   + f_4 * pc_x[k] * gpg_141[k];

        t_196[k] = f_4 * pc_z[k] * gpg_138[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, pc_x, pc_y, pc_z, fpg_50, fpg_145, \
                         gsg_55, gpf0_92, gpf1_92, gpg_140, gpg_141, \
                         gpg_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_10 * fpg_50[k]
                   + f_4 * pc_y[k] * gpg_140[k];

        t_198[k] = f_7 * gpf0_92[k]
                   - f_8 * gpf1_92[k]
                   + f_4 * pc_z[k] * gpg_140[k];

        t_199[k] = f_10 * fpg_145[k]
                   + f_1 * gsg_55[k]
                   + f_4 * pc_x[k] * gpg_145[k];

        t_200[k] = f_4 * pc_z[k] * gpg_141[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, pc_x, fpg_147, fpg_148, fpg_149, gsg_57, gsg_58, \
                         gsg_59, gpg_147, gpg_148, gpg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_10 * fpg_147[k]
                   + f_1 * gsg_57[k]
                   + f_4 * pc_x[k] * gpg_147[k];

        t_202[k] = f_10 * fpg_148[k]
                   + f_1 * gsg_58[k]
                   + f_4 * pc_x[k] * gpg_148[k];

        t_203[k] = f_10 * fpg_149[k]
                   + f_1 * gsg_59[k]
                   + f_4 * pc_x[k] * gpg_149[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pc_y, pc_z, fpg_55, gpf0_96, gpf0_97, \
                         gpf1_96, gpf1_97, gpg_145, gpg_146, gpg_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_10 * fpg_55[k]
                   + f_2 * gpf0_96[k]
                   - f_3 * gpf1_96[k]
                   + f_4 * pc_y[k] * gpg_145[k];

        t_205[k] = f_4 * pc_z[k] * gpg_145[k];

        t_206[k] = f_5 * gpf0_96[k]
                   - f_6 * gpf1_96[k]
                   + f_4 * pc_z[k] * gpg_146[k];

        t_207[k] = f_7 * gpf0_97[k]
                   - f_8 * gpf1_97[k]
                   + f_4 * pc_z[k] * gpg_147[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, pc_x, pc_y, pc_z, fpg_59, fpg_150, gpf0_99, \
                         gpf0_100, gpf1_99, gpf1_100, gpg_149, \
                         gpg_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_10 * fpg_59[k]
                   + f_4 * pc_y[k] * gpg_149[k];

        t_209[k] = f_2 * gpf0_99[k]
                   - f_3 * gpf1_99[k]
                   + f_4 * pc_z[k] * gpg_149[k];

        t_210[k] = f_10 * fpg_150[k]
                   + f_2 * gpf0_100[k]
                   - f_3 * gpf1_100[k]
                   + f_4 * pc_x[k] * gpg_150[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, pc_x, pc_y, pc_z, fpg_60, fpg_153, \
                         gsg_45, gpf0_103, gpf1_103, gpg_150, gpg_151, \
                         gpg_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = f_10 * fpg_60[k]
                   + f_1 * gsg_45[k]
                   + f_4 * pc_y[k] * gpg_150[k];

        t_212[k] = f_4 * pc_z[k] * gpg_150[k];

        t_213[k] = f_10 * fpg_153[k]
                   + f_7 * gpf0_103[k]
                   - f_8 * gpf1_103[k]
                   + f_4 * pc_x[k] * gpg_153[k];

        t_214[k] = f_4 * pc_z[k] * gpg_151[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, pc_x, pc_z, fpg_156, gpf0_100, gpf0_106, \
                         gpf1_100, gpf1_106, gpg_152, gpg_153, \
                         gpg_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_5 * gpf0_100[k]
                   - f_6 * gpf1_100[k]
                   + f_4 * pc_z[k] * gpg_152[k];

        t_216[k] = f_10 * fpg_156[k]
                   + f_5 * gpf0_106[k]
                   - f_6 * gpf1_106[k]
                   + f_4 * pc_x[k] * gpg_156[k];

        t_217[k] = f_4 * pc_z[k] * gpg_153[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, pc_x, pc_y, pc_z, fpg_65, fpg_160, \
                         gsg_50, gpf0_102, gpf1_102, gpg_155, gpg_156, \
                         gpg_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_10 * fpg_65[k]
                   + f_1 * gsg_50[k]
                   + f_4 * pc_y[k] * gpg_155[k];

        t_219[k] = f_7 * gpf0_102[k]
                   - f_8 * gpf1_102[k]
                   + f_4 * pc_z[k] * gpg_155[k];

        t_220[k] = f_10 * fpg_160[k]
                   + f_4 * pc_x[k] * gpg_160[k];

        t_221[k] = f_4 * pc_z[k] * gpg_156[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pa_x, pc_x, dph0_225, dph1_225, fph0_225, \
                         fpg_162, fpg_163, fpg_164, fph1_225, gpg_162, gpg_163, \
                         gpg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_10 * fpg_162[k]
                   + f_4 * pc_x[k] * gpg_162[k];

        t_223[k] = f_10 * fpg_163[k]
                   + f_4 * pc_x[k] * gpg_163[k];

        t_224[k] = f_10 * fpg_164[k]
                   + f_4 * pc_x[k] * gpg_164[k];

        t_225[k] = f_17 * dph0_225[k]
                   - f_18 * dph1_225[k]
                   + pa_x[k] * fph0_225[k]
                   - f_9 * pc_x[k] * fph1_225[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, pc_y, pc_z, fpg_74, gsg_59, gpf0_106, \
                         gpf0_107, gpf1_106, gpf1_107, gpg_160, gpg_161, gpg_162, \
                         gpg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_4 * pc_z[k] * gpg_160[k];

        t_227[k] = f_5 * gpf0_106[k]
                   - f_6 * gpf1_106[k]
                   + f_4 * pc_z[k] * gpg_161[k];

        t_228[k] = f_7 * gpf0_107[k]
                   - f_8 * gpf1_107[k]
                   + f_4 * pc_z[k] * gpg_162[k];

        t_229[k] = f_10 * fpg_74[k]
                   + f_1 * gsg_59[k]
                   + f_4 * pc_y[k] * gpg_164[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pb_z, pc_y, pc_z, fpg_75, gsh0_63, \
                         gsg_45, gsh1_63, gpf0_109, gpf1_109, gpg_164, \
                         gpg_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_2 * gpf0_109[k]
                   - f_3 * gpf1_109[k]
                   + f_4 * pc_z[k] * gpg_164[k];

        t_231[k] = pb_z[k] * gsh0_63[k]
                   - f_9 * pc_z[k] * gsh1_63[k];

        t_232[k] = f_10 * fpg_75[k]
                   + f_4 * pc_y[k] * gpg_165[k];

        t_233[k] = f_1 * gsg_45[k]
                   + f_4 * pc_z[k] * gpg_165[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, pb_z, pc_z, gsh0_66, gsh0_68, gsh0_69, \
                         gsg_46, gsg_47, gsh1_66, gsh1_68, gsh1_69, \
                         gpg_166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = pb_z[k] * gsh0_66[k]
                   - f_9 * pc_z[k] * gsh1_66[k];

        t_235[k] = f_1 * gsg_46[k]
                   + f_4 * pc_z[k] * gpg_166[k];

        t_236[k] = pb_z[k] * gsh0_68[k]
                   + f_10 * gsg_47[k]
                   - f_9 * pc_z[k] * gsh1_68[k];

        t_237[k] = pb_z[k] * gsh0_69[k]
                   - f_9 * pc_z[k] * gsh1_69[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, pb_z, pc_y, pc_z, fpg_80, gsh0_72, gsg_48, \
                         gsg_50, gsh1_72, gpg_168, gpg_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_1 * gsg_48[k]
                   + f_4 * pc_z[k] * gpg_168[k];

        t_239[k] = f_10 * fpg_80[k]
                   + f_4 * pc_y[k] * gpg_170[k];

        t_240[k] = pb_z[k] * gsh0_72[k]
                   + f_11 * gsg_50[k]
                   - f_9 * pc_z[k] * gsh1_72[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, pc_x, pc_z, fpg_175, fpg_177, fpg_178, \
                         gsg_51, gpg_171, gpg_175, gpg_177, gpg_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = f_10 * fpg_175[k]
                   + f_4 * pc_x[k] * gpg_175[k];

        t_242[k] = f_1 * gsg_51[k]
                   + f_4 * pc_z[k] * gpg_171[k];

        t_243[k] = f_10 * fpg_177[k]
                   + f_4 * pc_x[k] * gpg_177[k];

        t_244[k] = f_10 * fpg_178[k]
                   + f_4 * pc_x[k] * gpg_178[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, pb_z, pc_x, pc_z, fpg_179, gsh0_78, \
                         gsh0_80, gsg_55, gsg_56, gsh1_78, gsh1_80, gpg_175, \
                         gpg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_10 * fpg_179[k]
                   + f_4 * pc_x[k] * gpg_179[k];

        t_246[k] = pb_z[k] * gsh0_78[k]
                   - f_9 * pc_z[k] * gsh1_78[k];

        t_247[k] = f_1 * gsg_55[k]
                   + f_4 * pc_z[k] * gpg_175[k];

        t_248[k] = pb_z[k] * gsh0_80[k]
                   + f_10 * gsg_56[k]
                   - f_9 * pc_z[k] * gsh1_80[k];
    }
}

static auto
compute_prim_gph_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t dph0, const size_t dph1,
                                                          const size_t fph0, const size_t fpg,
                                                          const size_t fph1, const size_t gsh0,
                                                          const size_t gsg, const size_t gsh1,
                                                          const size_t gpf0, const size_t gpf1,
                                                          const size_t gpg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 2.0 / gamma;
    const auto f_3 = 2.0 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = gamma / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 2.5 / q;
    const auto f_15 = 1.5 / gamma;
    const auto f_16 = 1.5 * p / (gamma * q);
    const auto f_17 = 0.5 / p;
    const auto f_18 = 0.5 * gamma / (p * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dph0_0 = buffer.data(dph0 + 0);
    const auto *dph0_290 = buffer.data(dph0 + 290);
    const auto *dph0_291 = buffer.data(dph0 + 291);

    const auto *dph1_0 = buffer.data(dph1 + 0);
    const auto *dph1_290 = buffer.data(dph1 + 290);
    const auto *dph1_291 = buffer.data(dph1 + 291);

    const auto *fph0_66 = buffer.data(fph0 + 66);
    const auto *fph0_69 = buffer.data(fph0 + 69);
    const auto *fph0_73 = buffer.data(fph0 + 73);
    const auto *fph0_78 = buffer.data(fph0 + 78);
    const auto *fph0_85 = buffer.data(fph0 + 85);
    const auto *fph0_87 = buffer.data(fph0 + 87);
    const auto *fph0_90 = buffer.data(fph0 + 90);
    const auto *fph0_99 = buffer.data(fph0 + 99);
    const auto *fph0_126 = buffer.data(fph0 + 126);
    const auto *fph0_131 = buffer.data(fph0 + 131);
    const auto *fph0_135 = buffer.data(fph0 + 135);
    const auto *fph0_140 = buffer.data(fph0 + 140);
    const auto *fph0_146 = buffer.data(fph0 + 146);
    const auto *fph0_168 = buffer.data(fph0 + 168);
    const auto *fph0_170 = buffer.data(fph0 + 170);
    const auto *fph0_173 = buffer.data(fph0 + 173);
    const auto *fph0_177 = buffer.data(fph0 + 177);
    const auto *fph0_188 = buffer.data(fph0 + 188);
    const auto *fph0_290 = buffer.data(fph0 + 290);
    const auto *fph0_291 = buffer.data(fph0 + 291);

    const auto *fpg_45 = buffer.data(fpg + 45);
    const auto *fpg_48 = buffer.data(fpg + 48);
    const auto *fpg_55 = buffer.data(fpg + 55);
    const auto *fpg_60 = buffer.data(fpg + 60);
    const auto *fpg_63 = buffer.data(fpg + 63);
    const auto *fpg_70 = buffer.data(fpg + 70);
    const auto *fpg_74 = buffer.data(fpg + 74);
    const auto *fpg_78 = buffer.data(fpg + 78);
    const auto *fpg_85 = buffer.data(fpg + 85);
    const auto *fpg_89 = buffer.data(fpg + 89);
    const auto *fpg_90 = buffer.data(fpg + 90);
    const auto *fpg_92 = buffer.data(fpg + 92);
    const auto *fpg_95 = buffer.data(fpg + 95);
    const auto *fpg_102 = buffer.data(fpg + 102);
    const auto *fpg_103 = buffer.data(fpg + 103);
    const auto *fpg_104 = buffer.data(fpg + 104);
    const auto *fpg_105 = buffer.data(fpg + 105);
    const auto *fpg_107 = buffer.data(fpg + 107);
    const auto *fpg_110 = buffer.data(fpg + 110);
    const auto *fpg_119 = buffer.data(fpg + 119);
    const auto *fpg_120 = buffer.data(fpg + 120);
    const auto *fpg_122 = buffer.data(fpg + 122);
    const auto *fpg_125 = buffer.data(fpg + 125);
    const auto *fpg_130 = buffer.data(fpg + 130);
    const auto *fpg_132 = buffer.data(fpg + 132);
    const auto *fpg_133 = buffer.data(fpg + 133);
    const auto *fpg_134 = buffer.data(fpg + 134);
    const auto *fpg_191 = buffer.data(fpg + 191);
    const auto *fpg_192 = buffer.data(fpg + 192);
    const auto *fpg_193 = buffer.data(fpg + 193);
    const auto *fpg_195 = buffer.data(fpg + 195);
    const auto *fpg_200 = buffer.data(fpg + 200);
    const auto *fpg_204 = buffer.data(fpg + 204);
    const auto *fpg_205 = buffer.data(fpg + 205);
    const auto *fpg_206 = buffer.data(fpg + 206);
    const auto *fpg_207 = buffer.data(fpg + 207);
    const auto *fpg_208 = buffer.data(fpg + 208);
    const auto *fpg_209 = buffer.data(fpg + 209);
    const auto *fpg_213 = buffer.data(fpg + 213);
    const auto *fpg_216 = buffer.data(fpg + 216);
    const auto *fpg_220 = buffer.data(fpg + 220);
    const auto *fpg_221 = buffer.data(fpg + 221);
    const auto *fpg_222 = buffer.data(fpg + 222);
    const auto *fpg_223 = buffer.data(fpg + 223);
    const auto *fpg_224 = buffer.data(fpg + 224);
    const auto *fpg_230 = buffer.data(fpg + 230);
    const auto *fpg_234 = buffer.data(fpg + 234);
    const auto *fpg_235 = buffer.data(fpg + 235);
    const auto *fpg_236 = buffer.data(fpg + 236);
    const auto *fpg_237 = buffer.data(fpg + 237);
    const auto *fpg_239 = buffer.data(fpg + 239);
    const auto *fpg_250 = buffer.data(fpg + 250);
    const auto *fpg_251 = buffer.data(fpg + 251);
    const auto *fpg_252 = buffer.data(fpg + 252);
    const auto *fpg_254 = buffer.data(fpg + 254);
    const auto *fpg_255 = buffer.data(fpg + 255);

    const auto *fph1_66 = buffer.data(fph1 + 66);
    const auto *fph1_69 = buffer.data(fph1 + 69);
    const auto *fph1_73 = buffer.data(fph1 + 73);
    const auto *fph1_78 = buffer.data(fph1 + 78);
    const auto *fph1_85 = buffer.data(fph1 + 85);
    const auto *fph1_87 = buffer.data(fph1 + 87);
    const auto *fph1_90 = buffer.data(fph1 + 90);
    const auto *fph1_99 = buffer.data(fph1 + 99);
    const auto *fph1_126 = buffer.data(fph1 + 126);
    const auto *fph1_131 = buffer.data(fph1 + 131);
    const auto *fph1_135 = buffer.data(fph1 + 135);
    const auto *fph1_140 = buffer.data(fph1 + 140);
    const auto *fph1_146 = buffer.data(fph1 + 146);
    const auto *fph1_168 = buffer.data(fph1 + 168);
    const auto *fph1_170 = buffer.data(fph1 + 170);
    const auto *fph1_173 = buffer.data(fph1 + 173);
    const auto *fph1_177 = buffer.data(fph1 + 177);
    const auto *fph1_188 = buffer.data(fph1 + 188);
    const auto *fph1_290 = buffer.data(fph1 + 290);
    const auto *fph1_291 = buffer.data(fph1 + 291);

    const auto *gsh0_81 = buffer.data(gsh0 + 81);
    const auto *gsh0_83 = buffer.data(gsh0 + 83);
    const auto *gsh0_105 = buffer.data(gsh0 + 105);
    const auto *gsh0_108 = buffer.data(gsh0 + 108);
    const auto *gsh0_110 = buffer.data(gsh0 + 110);
    const auto *gsh0_111 = buffer.data(gsh0 + 111);
    const auto *gsh0_112 = buffer.data(gsh0 + 112);
    const auto *gsh0_114 = buffer.data(gsh0 + 114);
    const auto *gsh0_120 = buffer.data(gsh0 + 120);
    const auto *gsh0_121 = buffer.data(gsh0 + 121);
    const auto *gsh0_122 = buffer.data(gsh0 + 122);
    const auto *gsh0_123 = buffer.data(gsh0 + 123);
    const auto *gsh0_125 = buffer.data(gsh0 + 125);

    const auto *gsg_57 = buffer.data(gsg + 57);
    const auto *gsg_59 = buffer.data(gsg + 59);
    const auto *gsg_62 = buffer.data(gsg + 62);
    const auto *gsg_63 = buffer.data(gsg + 63);
    const auto *gsg_65 = buffer.data(gsg + 65);
    const auto *gsg_70 = buffer.data(gsg + 70);
    const auto *gsg_71 = buffer.data(gsg + 71);
    const auto *gsg_72 = buffer.data(gsg + 72);
    const auto *gsg_73 = buffer.data(gsg + 73);
    const auto *gsg_74 = buffer.data(gsg + 74);
    const auto *gsg_75 = buffer.data(gsg + 75);
    const auto *gsg_76 = buffer.data(gsg + 76);
    const auto *gsg_77 = buffer.data(gsg + 77);
    const auto *gsg_78 = buffer.data(gsg + 78);
    const auto *gsg_79 = buffer.data(gsg + 79);
    const auto *gsg_80 = buffer.data(gsg + 80);
    const auto *gsg_84 = buffer.data(gsg + 84);
    const auto *gsg_85 = buffer.data(gsg + 85);
    const auto *gsg_86 = buffer.data(gsg + 86);
    const auto *gsg_87 = buffer.data(gsg + 87);
    const auto *gsg_88 = buffer.data(gsg + 88);
    const auto *gsg_89 = buffer.data(gsg + 89);

    const auto *gsh1_81 = buffer.data(gsh1 + 81);
    const auto *gsh1_83 = buffer.data(gsh1 + 83);
    const auto *gsh1_105 = buffer.data(gsh1 + 105);
    const auto *gsh1_108 = buffer.data(gsh1 + 108);
    const auto *gsh1_110 = buffer.data(gsh1 + 110);
    const auto *gsh1_111 = buffer.data(gsh1 + 111);
    const auto *gsh1_112 = buffer.data(gsh1 + 112);
    const auto *gsh1_114 = buffer.data(gsh1 + 114);
    const auto *gsh1_120 = buffer.data(gsh1 + 120);
    const auto *gsh1_121 = buffer.data(gsh1 + 121);
    const auto *gsh1_122 = buffer.data(gsh1 + 122);
    const auto *gsh1_123 = buffer.data(gsh1 + 123);
    const auto *gsh1_125 = buffer.data(gsh1 + 125);

    const auto *gpf0_128 = buffer.data(gpf0 + 128);
    const auto *gpf0_129 = buffer.data(gpf0 + 129);
    const auto *gpf0_130 = buffer.data(gpf0 + 130);
    const auto *gpf0_135 = buffer.data(gpf0 + 135);
    const auto *gpf0_139 = buffer.data(gpf0 + 139);
    const auto *gpf0_143 = buffer.data(gpf0 + 143);
    const auto *gpf0_146 = buffer.data(gpf0 + 146);
    const auto *gpf0_148 = buffer.data(gpf0 + 148);
    const auto *gpf0_149 = buffer.data(gpf0 + 149);
    const auto *gpf0_150 = buffer.data(gpf0 + 150);
    const auto *gpf0_151 = buffer.data(gpf0 + 151);
    const auto *gpf0_152 = buffer.data(gpf0 + 152);
    const auto *gpf0_155 = buffer.data(gpf0 + 155);
    const auto *gpf0_156 = buffer.data(gpf0 + 156);
    const auto *gpf0_157 = buffer.data(gpf0 + 157);
    const auto *gpf0_158 = buffer.data(gpf0 + 158);
    const auto *gpf0_159 = buffer.data(gpf0 + 159);
    const auto *gpf0_170 = buffer.data(gpf0 + 170);

    const auto *gpf1_128 = buffer.data(gpf1 + 128);
    const auto *gpf1_129 = buffer.data(gpf1 + 129);
    const auto *gpf1_130 = buffer.data(gpf1 + 130);
    const auto *gpf1_135 = buffer.data(gpf1 + 135);
    const auto *gpf1_139 = buffer.data(gpf1 + 139);
    const auto *gpf1_143 = buffer.data(gpf1 + 143);
    const auto *gpf1_146 = buffer.data(gpf1 + 146);
    const auto *gpf1_148 = buffer.data(gpf1 + 148);
    const auto *gpf1_149 = buffer.data(gpf1 + 149);
    const auto *gpf1_150 = buffer.data(gpf1 + 150);
    const auto *gpf1_151 = buffer.data(gpf1 + 151);
    const auto *gpf1_152 = buffer.data(gpf1 + 152);
    const auto *gpf1_155 = buffer.data(gpf1 + 155);
    const auto *gpf1_156 = buffer.data(gpf1 + 156);
    const auto *gpf1_157 = buffer.data(gpf1 + 157);
    const auto *gpf1_158 = buffer.data(gpf1 + 158);
    const auto *gpf1_159 = buffer.data(gpf1 + 159);
    const auto *gpf1_170 = buffer.data(gpf1 + 170);

    const auto *gpg_179 = buffer.data(gpg + 179);
    const auto *gpg_180 = buffer.data(gpg + 180);
    const auto *gpg_182 = buffer.data(gpg + 182);
    const auto *gpg_183 = buffer.data(gpg + 183);
    const auto *gpg_185 = buffer.data(gpg + 185);
    const auto *gpg_190 = buffer.data(gpg + 190);
    const auto *gpg_191 = buffer.data(gpg + 191);
    const auto *gpg_192 = buffer.data(gpg + 192);
    const auto *gpg_193 = buffer.data(gpg + 193);
    const auto *gpg_194 = buffer.data(gpg + 194);
    const auto *gpg_195 = buffer.data(gpg + 195);
    const auto *gpg_197 = buffer.data(gpg + 197);
    const auto *gpg_198 = buffer.data(gpg + 198);
    const auto *gpg_200 = buffer.data(gpg + 200);
    const auto *gpg_204 = buffer.data(gpg + 204);
    const auto *gpg_205 = buffer.data(gpg + 205);
    const auto *gpg_206 = buffer.data(gpg + 206);
    const auto *gpg_207 = buffer.data(gpg + 207);
    const auto *gpg_208 = buffer.data(gpg + 208);
    const auto *gpg_209 = buffer.data(gpg + 209);
    const auto *gpg_210 = buffer.data(gpg + 210);
    const auto *gpg_212 = buffer.data(gpg + 212);
    const auto *gpg_213 = buffer.data(gpg + 213);
    const auto *gpg_215 = buffer.data(gpg + 215);
    const auto *gpg_216 = buffer.data(gpg + 216);
    const auto *gpg_220 = buffer.data(gpg + 220);
    const auto *gpg_221 = buffer.data(gpg + 221);
    const auto *gpg_222 = buffer.data(gpg + 222);
    const auto *gpg_223 = buffer.data(gpg + 223);
    const auto *gpg_224 = buffer.data(gpg + 224);
    const auto *gpg_225 = buffer.data(gpg + 225);
    const auto *gpg_226 = buffer.data(gpg + 226);
    const auto *gpg_227 = buffer.data(gpg + 227);
    const auto *gpg_228 = buffer.data(gpg + 228);
    const auto *gpg_229 = buffer.data(gpg + 229);
    const auto *gpg_230 = buffer.data(gpg + 230);
    const auto *gpg_234 = buffer.data(gpg + 234);
    const auto *gpg_235 = buffer.data(gpg + 235);
    const auto *gpg_236 = buffer.data(gpg + 236);
    const auto *gpg_237 = buffer.data(gpg + 237);
    const auto *gpg_238 = buffer.data(gpg + 238);
    const auto *gpg_239 = buffer.data(gpg + 239);
    const auto *gpg_240 = buffer.data(gpg + 240);
    const auto *gpg_242 = buffer.data(gpg + 242);
    const auto *gpg_245 = buffer.data(gpg + 245);
    const auto *gpg_249 = buffer.data(gpg + 249);
    const auto *gpg_250 = buffer.data(gpg + 250);
    const auto *gpg_251 = buffer.data(gpg + 251);
    const auto *gpg_252 = buffer.data(gpg + 252);
    const auto *gpg_254 = buffer.data(gpg + 254);
    const auto *gpg_255 = buffer.data(gpg + 255);
    const auto *gpg_256 = buffer.data(gpg + 256);
    const auto *gpg_257 = buffer.data(gpg + 257);

#pragma omp simd aligned(t_249, t_250, t_251, pb_z, pc_y, pc_z, fpg_89, gsh0_81, gsh0_83, \
                         gsg_57, gsg_59, gsh1_81, gsh1_83, gpg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = pb_z[k] * gsh0_81[k]
                   + f_11 * gsg_57[k]
                   - f_9 * pc_z[k] * gsh1_81[k];

        t_250[k] = f_10 * fpg_89[k]
                   + f_4 * pc_y[k] * gpg_179[k];

        t_251[k] = pb_z[k] * gsh0_83[k]
                   + f_12 * gsg_59[k]
                   - f_9 * pc_z[k] * gsh1_83[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pa_y, pa_z, pc_y, pc_z, fph0_66, \
                         fph0_126, fpg_45, fpg_90, fph1_66, fph1_126, \
                         gpg_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = pa_y[k] * fph0_126[k]
                   - f_9 * pc_y[k] * fph1_126[k];

        t_253[k] = f_1 * fpg_90[k]
                   + f_4 * pc_y[k] * gpg_180[k];

        t_254[k] = f_1 * fpg_45[k]
                   + f_4 * pc_z[k] * gpg_180[k];

        t_255[k] = pa_z[k] * fph0_66[k]
                   - f_9 * pc_z[k] * fph1_66[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, pa_y, pa_z, pc_y, pc_z, fph0_69, \
                         fph0_131, fpg_48, fpg_92, fph1_69, fph1_131, gpg_182, \
                         gpg_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_1 * fpg_92[k]
                   + f_4 * pc_y[k] * gpg_182[k];

        t_257[k] = pa_y[k] * fph0_131[k]
                   - f_9 * pc_y[k] * fph1_131[k];

        t_258[k] = pa_z[k] * fph0_69[k]
                   - f_9 * pc_z[k] * fph1_69[k];

        t_259[k] = f_1 * fpg_48[k]
                   + f_4 * pc_z[k] * gpg_183[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, pa_y, pa_z, pc_y, pc_z, fph0_73, fph0_135, \
                         fpg_95, fph1_73, fph1_135, gpg_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_1 * fpg_95[k]
                   + f_4 * pc_y[k] * gpg_185[k];

        t_261[k] = pa_y[k] * fph0_135[k]
                   - f_9 * pc_y[k] * fph1_135[k];

        t_262[k] = pa_z[k] * fph0_73[k]
                   - f_9 * pc_z[k] * fph1_73[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, pc_x, fpg_191, fpg_192, fpg_193, gsg_71, gsg_72, \
                         gsg_73, gpg_191, gpg_192, gpg_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_10 * fpg_191[k]
                   + f_1 * gsg_71[k]
                   + f_4 * pc_x[k] * gpg_191[k];

        t_264[k] = f_10 * fpg_192[k]
                   + f_1 * gsg_72[k]
                   + f_4 * pc_x[k] * gpg_192[k];

        t_265[k] = f_10 * fpg_193[k]
                   + f_1 * gsg_73[k]
                   + f_4 * pc_x[k] * gpg_193[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, pa_y, pa_z, pc_y, pc_z, fph0_78, fph0_140, \
                         fpg_55, fph1_78, fph1_140, gpg_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = pa_y[k] * fph0_140[k]
                   - f_9 * pc_y[k] * fph1_140[k];

        t_267[k] = pa_z[k] * fph0_78[k]
                   - f_9 * pc_z[k] * fph1_78[k];

        t_268[k] = f_1 * fpg_55[k]
                   + f_4 * pc_z[k] * gpg_190[k];
    }

#pragma omp simd aligned(t_269, t_270, t_271, pc_y, fpg_102, fpg_103, fpg_104, gpf0_128, \
                         gpf0_129, gpf1_128, gpf1_129, gpg_192, gpg_193, \
                         gpg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = f_1 * fpg_102[k]
                   + f_7 * gpf0_128[k]
                   - f_8 * gpf1_128[k]
                   + f_4 * pc_y[k] * gpg_192[k];

        t_270[k] = f_1 * fpg_103[k]
                   + f_5 * gpf0_129[k]
                   - f_6 * gpf1_129[k]
                   + f_4 * pc_y[k] * gpg_193[k];

        t_271[k] = f_1 * fpg_104[k]
                   + f_4 * pc_y[k] * gpg_194[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, pa_y, pa_z, pc_x, pc_y, pc_z, fph0_85, fph0_146, \
                         fpg_195, fph1_85, fph1_146, gpf0_130, gpf1_130, \
                         gpg_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = pa_y[k] * fph0_146[k]
                   - f_9 * pc_y[k] * fph1_146[k];

        t_273[k] = f_10 * fpg_195[k]
                   + f_2 * gpf0_130[k]
                   - f_3 * gpf1_130[k]
                   + f_4 * pc_x[k] * gpg_195[k];

        t_274[k] = pa_z[k] * fph0_85[k]
                   - f_9 * pc_z[k] * fph1_85[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, pa_z, pc_y, pc_z, fph0_87, fpg_60, fpg_107, \
                         fph1_87, gsg_62, gpg_195, gpg_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_1 * fpg_60[k]
                   + f_4 * pc_z[k] * gpg_195[k];

        t_276[k] = pa_z[k] * fph0_87[k]
                   - f_9 * pc_z[k] * fph1_87[k];

        t_277[k] = f_1 * fpg_107[k]
                   + f_1 * gsg_62[k]
                   + f_4 * pc_y[k] * gpg_197[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, pa_z, pc_x, pc_z, fph0_90, fpg_63, fpg_200, \
                         fph1_90, gpf0_135, gpf1_135, gpg_198, \
                         gpg_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_10 * fpg_200[k]
                   + f_7 * gpf0_135[k]
                   - f_8 * gpf1_135[k]
                   + f_4 * pc_x[k] * gpg_200[k];

        t_279[k] = pa_z[k] * fph0_90[k]
                   - f_9 * pc_z[k] * fph1_90[k];

        t_280[k] = f_1 * fpg_63[k]
                   + f_4 * pc_z[k] * gpg_198[k];
    }

#pragma omp simd aligned(t_281, t_282, t_283, pc_x, pc_y, fpg_110, fpg_204, fpg_205, gsg_65, \
                         gpf0_139, gpf1_139, gpg_200, gpg_204, \
                         gpg_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_281[k] = f_1 * fpg_110[k]
                   + f_1 * gsg_65[k]
                   + f_4 * pc_y[k] * gpg_200[k];

        t_282[k] = f_10 * fpg_204[k]
                   + f_5 * gpf0_139[k]
                   - f_6 * gpf1_139[k]
                   + f_4 * pc_x[k] * gpg_204[k];

        t_283[k] = f_10 * fpg_205[k]
                   + f_4 * pc_x[k] * gpg_205[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, pc_x, fpg_206, fpg_207, fpg_208, fpg_209, \
                         gpg_206, gpg_207, gpg_208, gpg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_10 * fpg_206[k]
                   + f_4 * pc_x[k] * gpg_206[k];

        t_285[k] = f_10 * fpg_207[k]
                   + f_4 * pc_x[k] * gpg_207[k];

        t_286[k] = f_10 * fpg_208[k]
                   + f_4 * pc_x[k] * gpg_208[k];

        t_287[k] = f_10 * fpg_209[k]
                   + f_4 * pc_x[k] * gpg_209[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, pa_x, pa_z, pc_x, pc_z, dph0_290, dph1_290, \
                         fph0_99, fph0_290, fpg_70, fph1_99, fph1_290, \
                         gpg_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = pa_z[k] * fph0_99[k]
                   - f_9 * pc_z[k] * fph1_99[k];

        t_289[k] = f_1 * fpg_70[k]
                   + f_4 * pc_z[k] * gpg_205[k];

        t_290[k] = f_17 * dph0_290[k]
                   - f_18 * dph1_290[k]
                   + pa_x[k] * fph0_290[k]
                   - f_9 * pc_x[k] * fph1_290[k];
    }

#pragma omp simd aligned(t_291, t_292, pa_x, pc_x, pc_y, dph0_291, dph1_291, fph0_291, \
                         fpg_119, fph1_291, gsg_74, gpg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_17 * dph0_291[k]
                   - f_18 * dph1_291[k]
                   + pa_x[k] * fph0_291[k]
                   - f_9 * pc_x[k] * fph1_291[k];

        t_292[k] = f_1 * fpg_119[k]
                   + f_1 * gsg_74[k]
                   + f_4 * pc_y[k] * gpg_209[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, pa_y, pc_y, pc_z, fph0_168, fpg_74, fpg_120, \
                         fph1_168, gpf0_139, gpf1_139, gpg_209, \
                         gpg_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = f_1 * fpg_74[k]
                   + f_2 * gpf0_139[k]
                   - f_3 * gpf1_139[k]
                   + f_4 * pc_z[k] * gpg_209[k];

        t_294[k] = pa_y[k] * fph0_168[k]
                   - f_9 * pc_y[k] * fph1_168[k];

        t_295[k] = f_1 * fpg_120[k]
                   + f_4 * pc_y[k] * gpg_210[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, pa_y, pc_x, pc_y, fph0_170, fpg_122, fpg_213, \
                         fph1_170, gpf0_143, gpf1_143, gpg_212, \
                         gpg_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = pa_y[k] * fph0_170[k]
                   - f_9 * pc_y[k] * fph1_170[k];

        t_297[k] = f_10 * fpg_213[k]
                   + f_7 * gpf0_143[k]
                   - f_8 * gpf1_143[k]
                   + f_4 * pc_x[k] * gpg_213[k];

        t_298[k] = f_1 * fpg_122[k]
                   + f_4 * pc_y[k] * gpg_212[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, pa_y, pc_x, pc_y, pc_z, fph0_173, fpg_78, \
                         fpg_216, fph1_173, gsg_63, gpf0_146, gpf1_146, gpg_213, \
                         gpg_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = pa_y[k] * fph0_173[k]
                   - f_9 * pc_y[k] * fph1_173[k];

        t_300[k] = f_10 * fpg_216[k]
                   + f_5 * gpf0_146[k]
                   - f_6 * gpf1_146[k]
                   + f_4 * pc_x[k] * gpg_216[k];

        t_301[k] = f_1 * fpg_78[k]
                   + f_1 * gsg_63[k]
                   + f_4 * pc_z[k] * gpg_213[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, t_305, pa_y, pc_x, pc_y, fph0_177, fpg_125, \
                         fpg_220, fpg_221, fph1_177, gpg_215, gpg_220, \
                         gpg_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = f_1 * fpg_125[k]
                   + f_4 * pc_y[k] * gpg_215[k];

        t_303[k] = pa_y[k] * fph0_177[k]
                   - f_9 * pc_y[k] * fph1_177[k];

        t_304[k] = f_10 * fpg_220[k]
                   + f_4 * pc_x[k] * gpg_220[k];

        t_305[k] = f_10 * fpg_221[k]
                   + f_4 * pc_x[k] * gpg_221[k];
    }

#pragma omp simd aligned(t_306, t_307, t_308, t_309, pc_x, pc_y, fpg_130, fpg_222, fpg_223, \
                         fpg_224, gpf0_146, gpf1_146, gpg_220, gpg_222, gpg_223, \
                         gpg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_306[k] = f_10 * fpg_222[k]
                   + f_4 * pc_x[k] * gpg_222[k];

        t_307[k] = f_10 * fpg_223[k]
                   + f_4 * pc_x[k] * gpg_223[k];

        t_308[k] = f_10 * fpg_224[k]
                   + f_4 * pc_x[k] * gpg_224[k];

        t_309[k] = f_1 * fpg_130[k]
                   + f_2 * gpf0_146[k]
                   - f_3 * gpf1_146[k]
                   + f_4 * pc_y[k] * gpg_220[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, pc_y, pc_z, fpg_85, fpg_132, fpg_133, gsg_70, \
                         gpf0_148, gpf0_149, gpf1_148, gpf1_149, gpg_220, gpg_222, \
                         gpg_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = f_1 * fpg_85[k]
                   + f_1 * gsg_70[k]
                   + f_4 * pc_z[k] * gpg_220[k];

        t_311[k] = f_1 * fpg_132[k]
                   + f_7 * gpf0_148[k]
                   - f_8 * gpf1_148[k]
                   + f_4 * pc_y[k] * gpg_222[k];

        t_312[k] = f_1 * fpg_133[k]
                   + f_5 * gpf0_149[k]
                   - f_6 * gpf1_149[k]
                   + f_4 * pc_y[k] * gpg_223[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, pa_y, pa_z, pc_y, pc_z, dph0_0, dph1_0, \
                         fph0_126, fph0_188, fpg_134, fph1_126, fph1_188, \
                         gpg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = f_1 * fpg_134[k]
                   + f_4 * pc_y[k] * gpg_224[k];

        t_314[k] = pa_y[k] * fph0_188[k]
                   - f_9 * pc_y[k] * fph1_188[k];

        t_315[k] = f_17 * dph0_0[k]
                   - f_18 * dph1_0[k]
                   + pa_z[k] * fph0_126[k]
                   - f_9 * pc_z[k] * fph1_126[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, t_319, pc_y, pc_z, fpg_90, gpf0_150, gpf1_150, \
                         gpg_225, gpg_226, gpg_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = f_4 * pc_y[k] * gpg_225[k];

        t_317[k] = f_10 * fpg_90[k]
                   + f_4 * pc_z[k] * gpg_225[k];

        t_318[k] = f_5 * gpf0_150[k]
                   - f_6 * gpf1_150[k]
                   + f_4 * pc_y[k] * gpg_226[k];

        t_319[k] = f_4 * pc_y[k] * gpg_227[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, pc_x, pc_y, fpg_230, gsg_80, gpf0_151, gpf0_152, \
                         gpf0_155, gpf1_151, gpf1_152, gpf1_155, gpg_228, gpg_229, \
                         gpg_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_10 * fpg_230[k]
                   + f_1 * gsg_80[k]
                   + f_7 * gpf0_155[k]
                   - f_8 * gpf1_155[k]
                   + f_4 * pc_x[k] * gpg_230[k];

        t_321[k] = f_7 * gpf0_151[k]
                   - f_8 * gpf1_151[k]
                   + f_4 * pc_y[k] * gpg_228[k];

        t_322[k] = f_5 * gpf0_152[k]
                   - f_6 * gpf1_152[k]
                   + f_4 * pc_y[k] * gpg_229[k];
    }

#pragma omp simd aligned(t_323, t_324, t_325, pc_x, pc_y, fpg_234, fpg_235, gsg_84, gsg_85, \
                         gpf0_159, gpf1_159, gpg_230, gpg_234, \
                         gpg_235 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_323[k] = f_4 * pc_y[k] * gpg_230[k];

        t_324[k] = f_10 * fpg_234[k]
                   + f_1 * gsg_84[k]
                   + f_5 * gpf0_159[k]
                   - f_6 * gpf1_159[k]
                   + f_4 * pc_x[k] * gpg_234[k];

        t_325[k] = f_10 * fpg_235[k]
                   + f_1 * gsg_85[k]
                   + f_4 * pc_x[k] * gpg_235[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, pc_x, pc_y, fpg_236, fpg_237, fpg_239, \
                         gsg_86, gsg_87, gsg_89, gpg_234, gpg_236, gpg_237, \
                         gpg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_10 * fpg_236[k]
                   + f_1 * gsg_86[k]
                   + f_4 * pc_x[k] * gpg_236[k];

        t_327[k] = f_10 * fpg_237[k]
                   + f_1 * gsg_87[k]
                   + f_4 * pc_x[k] * gpg_237[k];

        t_328[k] = f_4 * pc_y[k] * gpg_234[k];

        t_329[k] = f_10 * fpg_239[k]
                   + f_1 * gsg_89[k]
                   + f_4 * pc_x[k] * gpg_239[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, pc_y, gpf0_156, gpf0_157, gpf0_158, gpf1_156, \
                         gpf1_157, gpf1_158, gpg_235, gpg_236, \
                         gpg_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = f_2 * gpf0_156[k]
                   - f_3 * gpf1_156[k]
                   + f_4 * pc_y[k] * gpg_235[k];

        t_331[k] = f_15 * gpf0_157[k]
                   - f_16 * gpf1_157[k]
                   + f_4 * pc_y[k] * gpg_236[k];

        t_332[k] = f_7 * gpf0_158[k]
                   - f_8 * gpf1_158[k]
                   + f_4 * pc_y[k] * gpg_237[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, t_336, pb_y, pc_y, pc_z, fpg_104, gsh0_105, \
                         gsh1_105, gpf0_159, gpf1_159, gpg_238, \
                         gpg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = f_5 * gpf0_159[k]
                   - f_6 * gpf1_159[k]
                   + f_4 * pc_y[k] * gpg_238[k];

        t_334[k] = f_4 * pc_y[k] * gpg_239[k];

        t_335[k] = f_10 * fpg_104[k]
                   + f_2 * gpf0_159[k]
                   - f_3 * gpf1_159[k]
                   + f_4 * pc_z[k] * gpg_239[k];

        t_336[k] = pb_y[k] * gsh0_105[k]
                   - f_9 * pc_y[k] * gsh1_105[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, t_340, pb_y, pc_y, pc_z, fpg_105, gsh0_108, \
                         gsg_75, gsg_76, gsg_77, gsh1_108, gpg_240, \
                         gpg_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = f_1 * gsg_75[k]
                   + f_4 * pc_y[k] * gpg_240[k];

        t_338[k] = f_10 * fpg_105[k]
                   + f_4 * pc_z[k] * gpg_240[k];

        t_339[k] = pb_y[k] * gsh0_108[k]
                   + f_10 * gsg_76[k]
                   - f_9 * pc_y[k] * gsh1_108[k];

        t_340[k] = f_1 * gsg_77[k]
                   + f_4 * pc_y[k] * gpg_242[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, t_344, pb_y, pc_y, gsh0_110, gsh0_111, gsh0_112, \
                         gsg_78, gsg_79, gsg_80, gsh1_110, gsh1_111, gsh1_112, \
                         gpg_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = pb_y[k] * gsh0_110[k]
                   - f_9 * pc_y[k] * gsh1_110[k];

        t_342[k] = pb_y[k] * gsh0_111[k]
                   + f_11 * gsg_78[k]
                   - f_9 * pc_y[k] * gsh1_111[k];

        t_343[k] = pb_y[k] * gsh0_112[k]
                   + f_10 * gsg_79[k]
                   - f_9 * pc_y[k] * gsh1_112[k];

        t_344[k] = f_1 * gsg_80[k]
                   + f_4 * pc_y[k] * gpg_245[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, pb_y, pc_x, pc_y, fpg_250, fpg_251, \
                         fpg_252, gsh0_114, gsh1_114, gpg_250, gpg_251, \
                         gpg_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = pb_y[k] * gsh0_114[k]
                   - f_9 * pc_y[k] * gsh1_114[k];

        t_346[k] = f_10 * fpg_250[k]
                   + f_4 * pc_x[k] * gpg_250[k];

        t_347[k] = f_10 * fpg_251[k]
                   + f_4 * pc_x[k] * gpg_251[k];

        t_348[k] = f_10 * fpg_252[k]
                   + f_4 * pc_x[k] * gpg_252[k];
    }

#pragma omp simd aligned(t_349, t_350, t_351, pb_y, pc_x, pc_y, fpg_254, gsh0_120, gsg_84, \
                         gsg_85, gsh1_120, gpg_249, gpg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_349[k] = f_1 * gsg_84[k]
                   + f_4 * pc_y[k] * gpg_249[k];

        t_350[k] = f_10 * fpg_254[k]
                   + f_4 * pc_x[k] * gpg_254[k];

        t_351[k] = pb_y[k] * gsh0_120[k]
                   + f_12 * gsg_85[k]
                   - f_9 * pc_y[k] * gsh1_120[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, pb_y, pc_y, gsh0_121, gsh0_122, gsh0_123, \
                         gsg_86, gsg_87, gsg_88, gsh1_121, gsh1_122, \
                         gsh1_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = pb_y[k] * gsh0_121[k]
                   + f_0 * gsg_86[k]
                   - f_9 * pc_y[k] * gsh1_121[k];

        t_353[k] = pb_y[k] * gsh0_122[k]
                   + f_11 * gsg_87[k]
                   - f_9 * pc_y[k] * gsh1_122[k];

        t_354[k] = pb_y[k] * gsh0_123[k]
                   + f_10 * gsg_88[k]
                   - f_9 * pc_y[k] * gsh1_123[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, pb_y, pc_x, pc_y, fpg_255, gsh0_125, \
                         gsg_89, gsh1_125, gpf0_170, gpf1_170, gpg_254, \
                         gpg_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = f_1 * gsg_89[k]
                   + f_4 * pc_y[k] * gpg_254[k];

        t_356[k] = pb_y[k] * gsh0_125[k]
                   - f_9 * pc_y[k] * gsh1_125[k];

        t_357[k] = f_10 * fpg_255[k]
                   + f_2 * gpf0_170[k]
                   - f_3 * gpf1_170[k]
                   + f_4 * pc_x[k] * gpg_255[k];

        t_358[k] = f_4 * pc_y[k] * gpg_255[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, pc_y, pc_z, fpg_120, gsg_75, gpf0_170, gpf1_170, \
                         gpg_255, gpg_256, gpg_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_10 * fpg_120[k]
                   + f_1 * gsg_75[k]
                   + f_4 * pc_z[k] * gpg_255[k];

        t_360[k] = f_5 * gpf0_170[k]
                   - f_6 * gpf1_170[k]
                   + f_4 * pc_y[k] * gpg_256[k];

        t_361[k] = f_4 * pc_y[k] * gpg_257[k];
    }
}

static auto
compute_prim_gph_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t dph0, const size_t dph1,
                                                          const size_t fph0, const size_t fpg,
                                                          const size_t fph1, const size_t gsh0,
                                                          const size_t gsg, const size_t gsh1,
                                                          const size_t gpf0, const size_t gpf1,
                                                          const size_t gpg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 0.5 / q;
    const auto f_2 = 2.0 / gamma;
    const auto f_3 = 2.0 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = gamma / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 2.5 / q;
    const auto f_15 = 1.5 / gamma;
    const auto f_16 = 1.5 * p / (gamma * q);
    const auto f_17 = 0.5 / p;
    const auto f_18 = 0.5 * gamma / (p * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dph0_377 = buffer.data(dph0 + 377);

    const auto *dph1_377 = buffer.data(dph1 + 377);

    const auto *fph0_189 = buffer.data(fph0 + 189);
    const auto *fph0_192 = buffer.data(fph0 + 192);
    const auto *fph0_195 = buffer.data(fph0 + 195);
    const auto *fph0_199 = buffer.data(fph0 + 199);
    const auto *fph0_210 = buffer.data(fph0 + 210);
    const auto *fph0_211 = buffer.data(fph0 + 211);
    const auto *fph0_213 = buffer.data(fph0 + 213);
    const auto *fph0_216 = buffer.data(fph0 + 216);
    const auto *fph0_377 = buffer.data(fph0 + 377);
    const auto *fph0_399 = buffer.data(fph0 + 399);
    const auto *fph0_402 = buffer.data(fph0 + 402);
    const auto *fph0_405 = buffer.data(fph0 + 405);
    const auto *fph0_414 = buffer.data(fph0 + 414);
    const auto *fph0_416 = buffer.data(fph0 + 416);
    const auto *fph0_417 = buffer.data(fph0 + 417);
    const auto *fph0_418 = buffer.data(fph0 + 418);
    const auto *fph0_419 = buffer.data(fph0 + 419);
    const auto *fph0_425 = buffer.data(fph0 + 425);
    const auto *fph0_429 = buffer.data(fph0 + 429);
    const auto *fph0_435 = buffer.data(fph0 + 435);
    const auto *fph0_437 = buffer.data(fph0 + 437);
    const auto *fph0_438 = buffer.data(fph0 + 438);
    const auto *fph0_440 = buffer.data(fph0 + 440);
    const auto *fph0_467 = buffer.data(fph0 + 467);
    const auto *fph0_471 = buffer.data(fph0 + 471);

    const auto *fpg_135 = buffer.data(fpg + 135);
    const auto *fpg_138 = buffer.data(fpg + 138);
    const auto *fpg_140 = buffer.data(fpg + 140);
    const auto *fpg_145 = buffer.data(fpg + 145);
    const auto *fpg_149 = buffer.data(fpg + 149);
    const auto *fpg_150 = buffer.data(fpg + 150);
    const auto *fpg_153 = buffer.data(fpg + 153);
    const auto *fpg_155 = buffer.data(fpg + 155);
    const auto *fpg_165 = buffer.data(fpg + 165);
    const auto *fpg_170 = buffer.data(fpg + 170);
    const auto *fpg_179 = buffer.data(fpg + 179);
    const auto *fpg_180 = buffer.data(fpg + 180);
    const auto *fpg_182 = buffer.data(fpg + 182);
    const auto *fpg_185 = buffer.data(fpg + 185);
    const auto *fpg_190 = buffer.data(fpg + 190);
    const auto *fpg_192 = buffer.data(fpg + 192);
    const auto *fpg_193 = buffer.data(fpg + 193);
    const auto *fpg_194 = buffer.data(fpg + 194);
    const auto *fpg_197 = buffer.data(fpg + 197);
    const auto *fpg_200 = buffer.data(fpg + 200);
    const auto *fpg_260 = buffer.data(fpg + 260);
    const auto *fpg_264 = buffer.data(fpg + 264);
    const auto *fpg_265 = buffer.data(fpg + 265);
    const auto *fpg_266 = buffer.data(fpg + 266);
    const auto *fpg_267 = buffer.data(fpg + 267);
    const auto *fpg_269 = buffer.data(fpg + 269);
    const auto *fpg_270 = buffer.data(fpg + 270);
    const auto *fpg_273 = buffer.data(fpg + 273);
    const auto *fpg_276 = buffer.data(fpg + 276);
    const auto *fpg_280 = buffer.data(fpg + 280);
    const auto *fpg_282 = buffer.data(fpg + 282);
    const auto *fpg_283 = buffer.data(fpg + 283);
    const auto *fpg_284 = buffer.data(fpg + 284);
    const auto *fpg_285 = buffer.data(fpg + 285);
    const auto *fpg_288 = buffer.data(fpg + 288);
    const auto *fpg_291 = buffer.data(fpg + 291);
    const auto *fpg_295 = buffer.data(fpg + 295);
    const auto *fpg_297 = buffer.data(fpg + 297);
    const auto *fpg_298 = buffer.data(fpg + 298);
    const auto *fpg_299 = buffer.data(fpg + 299);
    const auto *fpg_305 = buffer.data(fpg + 305);
    const auto *fpg_309 = buffer.data(fpg + 309);
    const auto *fpg_310 = buffer.data(fpg + 310);
    const auto *fpg_312 = buffer.data(fpg + 312);
    const auto *fpg_313 = buffer.data(fpg + 313);
    const auto *fpg_314 = buffer.data(fpg + 314);
    const auto *fpg_320 = buffer.data(fpg + 320);
    const auto *fpg_324 = buffer.data(fpg + 324);
    const auto *fpg_326 = buffer.data(fpg + 326);
    const auto *fpg_327 = buffer.data(fpg + 327);
    const auto *fpg_328 = buffer.data(fpg + 328);
    const auto *fpg_329 = buffer.data(fpg + 329);
    const auto *fpg_335 = buffer.data(fpg + 335);
    const auto *fpg_339 = buffer.data(fpg + 339);
    const auto *fpg_340 = buffer.data(fpg + 340);
    const auto *fpg_341 = buffer.data(fpg + 341);
    const auto *fpg_342 = buffer.data(fpg + 342);
    const auto *fpg_343 = buffer.data(fpg + 343);
    const auto *fpg_344 = buffer.data(fpg + 344);

    const auto *fph1_189 = buffer.data(fph1 + 189);
    const auto *fph1_192 = buffer.data(fph1 + 192);
    const auto *fph1_195 = buffer.data(fph1 + 195);
    const auto *fph1_199 = buffer.data(fph1 + 199);
    const auto *fph1_210 = buffer.data(fph1 + 210);
    const auto *fph1_211 = buffer.data(fph1 + 211);
    const auto *fph1_213 = buffer.data(fph1 + 213);
    const auto *fph1_216 = buffer.data(fph1 + 216);
    const auto *fph1_377 = buffer.data(fph1 + 377);
    const auto *fph1_399 = buffer.data(fph1 + 399);
    const auto *fph1_402 = buffer.data(fph1 + 402);
    const auto *fph1_405 = buffer.data(fph1 + 405);
    const auto *fph1_414 = buffer.data(fph1 + 414);
    const auto *fph1_416 = buffer.data(fph1 + 416);
    const auto *fph1_417 = buffer.data(fph1 + 417);
    const auto *fph1_418 = buffer.data(fph1 + 418);
    const auto *fph1_419 = buffer.data(fph1 + 419);
    const auto *fph1_425 = buffer.data(fph1 + 425);
    const auto *fph1_429 = buffer.data(fph1 + 429);
    const auto *fph1_435 = buffer.data(fph1 + 435);
    const auto *fph1_437 = buffer.data(fph1 + 437);
    const auto *fph1_438 = buffer.data(fph1 + 438);
    const auto *fph1_440 = buffer.data(fph1 + 440);
    const auto *fph1_467 = buffer.data(fph1 + 467);
    const auto *fph1_471 = buffer.data(fph1 + 471);

    const auto *gsh0_126 = buffer.data(gsh0 + 126);
    const auto *gsh0_129 = buffer.data(gsh0 + 129);
    const auto *gsh0_132 = buffer.data(gsh0 + 132);

    const auto *gsg_90 = buffer.data(gsg + 90);
    const auto *gsg_91 = buffer.data(gsg + 91);
    const auto *gsg_93 = buffer.data(gsg + 93);
    const auto *gsg_95 = buffer.data(gsg + 95);
    const auto *gsg_96 = buffer.data(gsg + 96);
    const auto *gsg_100 = buffer.data(gsg + 100);
    const auto *gsg_102 = buffer.data(gsg + 102);
    const auto *gsg_103 = buffer.data(gsg + 103);
    const auto *gsg_104 = buffer.data(gsg + 104);
    const auto *gsg_107 = buffer.data(gsg + 107);
    const auto *gsg_110 = buffer.data(gsg + 110);
    const auto *gsg_114 = buffer.data(gsg + 114);
    const auto *gsg_116 = buffer.data(gsg + 116);
    const auto *gsg_117 = buffer.data(gsg + 117);
    const auto *gsg_118 = buffer.data(gsg + 118);
    const auto *gsg_119 = buffer.data(gsg + 119);

    const auto *gsh1_126 = buffer.data(gsh1 + 126);
    const auto *gsh1_129 = buffer.data(gsh1 + 129);
    const auto *gsh1_132 = buffer.data(gsh1 + 132);

    const auto *gpf0_171 = buffer.data(gpf0 + 171);
    const auto *gpf0_172 = buffer.data(gpf0 + 172);
    const auto *gpf0_175 = buffer.data(gpf0 + 175);
    const auto *gpf0_176 = buffer.data(gpf0 + 176);
    const auto *gpf0_177 = buffer.data(gpf0 + 177);
    const auto *gpf0_178 = buffer.data(gpf0 + 178);
    const auto *gpf0_179 = buffer.data(gpf0 + 179);
    const auto *gpf0_180 = buffer.data(gpf0 + 180);
    const auto *gpf0_182 = buffer.data(gpf0 + 182);
    const auto *gpf0_183 = buffer.data(gpf0 + 183);
    const auto *gpf0_186 = buffer.data(gpf0 + 186);
    const auto *gpf0_187 = buffer.data(gpf0 + 187);
    const auto *gpf0_189 = buffer.data(gpf0 + 189);
    const auto *gpf0_190 = buffer.data(gpf0 + 190);
    const auto *gpf0_192 = buffer.data(gpf0 + 192);
    const auto *gpf0_215 = buffer.data(gpf0 + 215);
    const auto *gpf0_216 = buffer.data(gpf0 + 216);
    const auto *gpf0_218 = buffer.data(gpf0 + 218);
    const auto *gpf0_219 = buffer.data(gpf0 + 219);

    const auto *gpf1_171 = buffer.data(gpf1 + 171);
    const auto *gpf1_172 = buffer.data(gpf1 + 172);
    const auto *gpf1_175 = buffer.data(gpf1 + 175);
    const auto *gpf1_176 = buffer.data(gpf1 + 176);
    const auto *gpf1_177 = buffer.data(gpf1 + 177);
    const auto *gpf1_178 = buffer.data(gpf1 + 178);
    const auto *gpf1_179 = buffer.data(gpf1 + 179);
    const auto *gpf1_180 = buffer.data(gpf1 + 180);
    const auto *gpf1_182 = buffer.data(gpf1 + 182);
    const auto *gpf1_183 = buffer.data(gpf1 + 183);
    const auto *gpf1_186 = buffer.data(gpf1 + 186);
    const auto *gpf1_187 = buffer.data(gpf1 + 187);
    const auto *gpf1_189 = buffer.data(gpf1 + 189);
    const auto *gpf1_190 = buffer.data(gpf1 + 190);
    const auto *gpf1_192 = buffer.data(gpf1 + 192);
    const auto *gpf1_215 = buffer.data(gpf1 + 215);
    const auto *gpf1_216 = buffer.data(gpf1 + 216);
    const auto *gpf1_218 = buffer.data(gpf1 + 218);
    const auto *gpf1_219 = buffer.data(gpf1 + 219);

    const auto *gpg_258 = buffer.data(gpg + 258);
    const auto *gpg_259 = buffer.data(gpg + 259);
    const auto *gpg_260 = buffer.data(gpg + 260);
    const auto *gpg_264 = buffer.data(gpg + 264);
    const auto *gpg_265 = buffer.data(gpg + 265);
    const auto *gpg_266 = buffer.data(gpg + 266);
    const auto *gpg_267 = buffer.data(gpg + 267);
    const auto *gpg_268 = buffer.data(gpg + 268);
    const auto *gpg_269 = buffer.data(gpg + 269);
    const auto *gpg_270 = buffer.data(gpg + 270);
    const auto *gpg_271 = buffer.data(gpg + 271);
    const auto *gpg_272 = buffer.data(gpg + 272);
    const auto *gpg_273 = buffer.data(gpg + 273);
    const auto *gpg_275 = buffer.data(gpg + 275);
    const auto *gpg_276 = buffer.data(gpg + 276);
    const auto *gpg_280 = buffer.data(gpg + 280);
    const auto *gpg_281 = buffer.data(gpg + 281);
    const auto *gpg_282 = buffer.data(gpg + 282);
    const auto *gpg_283 = buffer.data(gpg + 283);
    const auto *gpg_284 = buffer.data(gpg + 284);
    const auto *gpg_285 = buffer.data(gpg + 285);
    const auto *gpg_286 = buffer.data(gpg + 286);
    const auto *gpg_287 = buffer.data(gpg + 287);
    const auto *gpg_288 = buffer.data(gpg + 288);
    const auto *gpg_290 = buffer.data(gpg + 290);
    const auto *gpg_291 = buffer.data(gpg + 291);
    const auto *gpg_295 = buffer.data(gpg + 295);
    const auto *gpg_297 = buffer.data(gpg + 297);
    const auto *gpg_298 = buffer.data(gpg + 298);
    const auto *gpg_299 = buffer.data(gpg + 299);
    const auto *gpg_300 = buffer.data(gpg + 300);
    const auto *gpg_301 = buffer.data(gpg + 301);
    const auto *gpg_303 = buffer.data(gpg + 303);
    const auto *gpg_305 = buffer.data(gpg + 305);
    const auto *gpg_306 = buffer.data(gpg + 306);
    const auto *gpg_310 = buffer.data(gpg + 310);
    const auto *gpg_312 = buffer.data(gpg + 312);
    const auto *gpg_313 = buffer.data(gpg + 313);
    const auto *gpg_314 = buffer.data(gpg + 314);
    const auto *gpg_315 = buffer.data(gpg + 315);
    const auto *gpg_317 = buffer.data(gpg + 317);
    const auto *gpg_318 = buffer.data(gpg + 318);
    const auto *gpg_320 = buffer.data(gpg + 320);
    const auto *gpg_324 = buffer.data(gpg + 324);
    const auto *gpg_325 = buffer.data(gpg + 325);
    const auto *gpg_326 = buffer.data(gpg + 326);
    const auto *gpg_327 = buffer.data(gpg + 327);
    const auto *gpg_328 = buffer.data(gpg + 328);
    const auto *gpg_329 = buffer.data(gpg + 329);
    const auto *gpg_330 = buffer.data(gpg + 330);
    const auto *gpg_332 = buffer.data(gpg + 332);
    const auto *gpg_333 = buffer.data(gpg + 333);
    const auto *gpg_335 = buffer.data(gpg + 335);
    const auto *gpg_340 = buffer.data(gpg + 340);
    const auto *gpg_341 = buffer.data(gpg + 341);
    const auto *gpg_342 = buffer.data(gpg + 342);
    const auto *gpg_343 = buffer.data(gpg + 343);
    const auto *gpg_344 = buffer.data(gpg + 344);

#pragma omp simd aligned(t_362, t_363, t_364, t_365, pc_x, pc_y, fpg_260, gpf0_171, gpf0_172, \
                         gpf0_175, gpf1_171, gpf1_172, gpf1_175, gpg_258, gpg_259, \
                         gpg_260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_10 * fpg_260[k]
                   + f_7 * gpf0_175[k]
                   - f_8 * gpf1_175[k]
                   + f_4 * pc_x[k] * gpg_260[k];

        t_363[k] = f_7 * gpf0_171[k]
                   - f_8 * gpf1_171[k]
                   + f_4 * pc_y[k] * gpg_258[k];

        t_364[k] = f_5 * gpf0_172[k]
                   - f_6 * gpf1_172[k]
                   + f_4 * pc_y[k] * gpg_259[k];

        t_365[k] = f_4 * pc_y[k] * gpg_260[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pc_x, fpg_264, fpg_265, fpg_266, fpg_267, \
                         gpf0_179, gpf1_179, gpg_264, gpg_265, gpg_266, \
                         gpg_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_10 * fpg_264[k]
                   + f_5 * gpf0_179[k]
                   - f_6 * gpf1_179[k]
                   + f_4 * pc_x[k] * gpg_264[k];

        t_367[k] = f_10 * fpg_265[k]
                   + f_4 * pc_x[k] * gpg_265[k];

        t_368[k] = f_10 * fpg_266[k]
                   + f_4 * pc_x[k] * gpg_266[k];

        t_369[k] = f_10 * fpg_267[k]
                   + f_4 * pc_x[k] * gpg_267[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, pc_x, pc_y, fpg_269, gpf0_176, gpf0_177, \
                         gpf1_176, gpf1_177, gpg_264, gpg_265, gpg_266, \
                         gpg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = f_4 * pc_y[k] * gpg_264[k];

        t_371[k] = f_10 * fpg_269[k]
                   + f_4 * pc_x[k] * gpg_269[k];

        t_372[k] = f_2 * gpf0_176[k]
                   - f_3 * gpf1_176[k]
                   + f_4 * pc_y[k] * gpg_265[k];

        t_373[k] = f_15 * gpf0_177[k]
                   - f_16 * gpf1_177[k]
                   + f_4 * pc_y[k] * gpg_266[k];
    }

#pragma omp simd aligned(t_374, t_375, t_376, pc_y, gpf0_178, gpf0_179, gpf1_178, gpf1_179, \
                         gpg_267, gpg_268, gpg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = f_7 * gpf0_178[k]
                   - f_8 * gpf1_178[k]
                   + f_4 * pc_y[k] * gpg_267[k];

        t_375[k] = f_5 * gpf0_179[k]
                   - f_6 * gpf1_179[k]
                   + f_4 * pc_y[k] * gpg_268[k];

        t_376[k] = f_4 * pc_y[k] * gpg_269[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, pa_x, pc_x, pc_y, dph0_377, dph1_377, fph0_377, \
                         fpg_135, fpg_270, fph1_377, gsg_90, gpf0_180, gpf1_180, \
                         gpg_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = f_17 * dph0_377[k]
                   - f_18 * dph1_377[k]
                   + pa_x[k] * fph0_377[k]
                   - f_9 * pc_x[k] * fph1_377[k];

        t_378[k] = f_1 * fpg_270[k]
                   + f_1 * gsg_90[k]
                   + f_2 * gpf0_180[k]
                   - f_3 * gpf1_180[k]
                   + f_4 * pc_x[k] * gpg_270[k];

        t_379[k] = f_11 * fpg_135[k]
                   + f_4 * pc_y[k] * gpg_270[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, pc_x, pc_z, fpg_273, gsg_93, gpf0_180, \
                         gpf0_183, gpf1_180, gpf1_183, gpg_270, gpg_271, gpg_272, \
                         gpg_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = f_4 * pc_z[k] * gpg_270[k];

        t_381[k] = f_1 * fpg_273[k]
                   + f_1 * gsg_93[k]
                   + f_7 * gpf0_183[k]
                   - f_8 * gpf1_183[k]
                   + f_4 * pc_x[k] * gpg_273[k];

        t_382[k] = f_4 * pc_z[k] * gpg_271[k];

        t_383[k] = f_5 * gpf0_180[k]
                   - f_6 * gpf1_180[k]
                   + f_4 * pc_z[k] * gpg_272[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, pc_x, pc_y, pc_z, fpg_140, fpg_276, gsg_96, \
                         gpf0_186, gpf1_186, gpg_273, gpg_275, \
                         gpg_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = f_1 * fpg_276[k]
                   + f_1 * gsg_96[k]
                   + f_5 * gpf0_186[k]
                   - f_6 * gpf1_186[k]
                   + f_4 * pc_x[k] * gpg_276[k];

        t_385[k] = f_4 * pc_z[k] * gpg_273[k];

        t_386[k] = f_11 * fpg_140[k]
                   + f_4 * pc_y[k] * gpg_275[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, pc_x, pc_z, fpg_280, fpg_282, gsg_100, \
                         gsg_102, gpf0_182, gpf1_182, gpg_275, gpg_276, gpg_280, \
                         gpg_282 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = f_7 * gpf0_182[k]
                   - f_8 * gpf1_182[k]
                   + f_4 * pc_z[k] * gpg_275[k];

        t_388[k] = f_1 * fpg_280[k]
                   + f_1 * gsg_100[k]
                   + f_4 * pc_x[k] * gpg_280[k];

        t_389[k] = f_4 * pc_z[k] * gpg_276[k];

        t_390[k] = f_1 * fpg_282[k]
                   + f_1 * gsg_102[k]
                   + f_4 * pc_x[k] * gpg_282[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, pc_x, pc_y, fpg_145, fpg_283, fpg_284, gsg_103, \
                         gsg_104, gpf0_186, gpf1_186, gpg_280, gpg_283, \
                         gpg_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_1 * fpg_283[k]
                   + f_1 * gsg_103[k]
                   + f_4 * pc_x[k] * gpg_283[k];

        t_392[k] = f_1 * fpg_284[k]
                   + f_1 * gsg_104[k]
                   + f_4 * pc_x[k] * gpg_284[k];

        t_393[k] = f_11 * fpg_145[k]
                   + f_2 * gpf0_186[k]
                   - f_3 * gpf1_186[k]
                   + f_4 * pc_y[k] * gpg_280[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, t_397, pc_y, pc_z, fpg_149, gpf0_186, gpf0_187, \
                         gpf1_186, gpf1_187, gpg_280, gpg_281, gpg_282, \
                         gpg_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = f_4 * pc_z[k] * gpg_280[k];

        t_395[k] = f_5 * gpf0_186[k]
                   - f_6 * gpf1_186[k]
                   + f_4 * pc_z[k] * gpg_281[k];

        t_396[k] = f_7 * gpf0_187[k]
                   - f_8 * gpf1_187[k]
                   + f_4 * pc_z[k] * gpg_282[k];

        t_397[k] = f_11 * fpg_149[k]
                   + f_4 * pc_y[k] * gpg_284[k];
    }

#pragma omp simd aligned(t_398, t_399, t_400, pa_x, pc_x, pc_y, pc_z, fph0_399, fpg_150, \
                         fpg_285, fph1_399, gsg_90, gpf0_189, gpf1_189, gpg_284, \
                         gpg_285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_398[k] = f_2 * gpf0_189[k]
                   - f_3 * gpf1_189[k]
                   + f_4 * pc_z[k] * gpg_284[k];

        t_399[k] = pa_x[k] * fph0_399[k]
                   + f_12 * fpg_285[k]
                   - f_9 * pc_x[k] * fph1_399[k];

        t_400[k] = f_11 * fpg_150[k]
                   + f_1 * gsg_90[k]
                   + f_4 * pc_y[k] * gpg_285[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, t_404, pa_x, pc_x, pc_z, fph0_402, fpg_288, \
                         fph1_402, gpf0_190, gpf1_190, gpg_285, gpg_286, \
                         gpg_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_4 * pc_z[k] * gpg_285[k];

        t_402[k] = pa_x[k] * fph0_402[k]
                   + f_11 * fpg_288[k]
                   - f_9 * pc_x[k] * fph1_402[k];

        t_403[k] = f_4 * pc_z[k] * gpg_286[k];

        t_404[k] = f_5 * gpf0_190[k]
                   - f_6 * gpf1_190[k]
                   + f_4 * pc_z[k] * gpg_287[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, pa_x, pc_x, pc_y, pc_z, fph0_405, fpg_155, \
                         fpg_291, fph1_405, gsg_95, gpg_288, gpg_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = pa_x[k] * fph0_405[k]
                   + f_10 * fpg_291[k]
                   - f_9 * pc_x[k] * fph1_405[k];

        t_406[k] = f_4 * pc_z[k] * gpg_288[k];

        t_407[k] = f_11 * fpg_155[k]
                   + f_1 * gsg_95[k]
                   + f_4 * pc_y[k] * gpg_290[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, t_411, pc_x, pc_z, fpg_295, fpg_297, gpf0_192, \
                         gpf1_192, gpg_290, gpg_291, gpg_295, gpg_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = f_7 * gpf0_192[k]
                   - f_8 * gpf1_192[k]
                   + f_4 * pc_z[k] * gpg_290[k];

        t_409[k] = f_1 * fpg_295[k]
                   + f_4 * pc_x[k] * gpg_295[k];

        t_410[k] = f_4 * pc_z[k] * gpg_291[k];

        t_411[k] = f_1 * fpg_297[k]
                   + f_4 * pc_x[k] * gpg_297[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, t_415, pa_x, pc_x, pc_z, fph0_414, fpg_298, \
                         fpg_299, fph1_414, gpg_295, gpg_298, gpg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = f_1 * fpg_298[k]
                   + f_4 * pc_x[k] * gpg_298[k];

        t_413[k] = f_1 * fpg_299[k]
                   + f_4 * pc_x[k] * gpg_299[k];

        t_414[k] = pa_x[k] * fph0_414[k]
                   - f_9 * pc_x[k] * fph1_414[k];

        t_415[k] = f_4 * pc_z[k] * gpg_295[k];
    }

#pragma omp simd aligned(t_416, t_417, t_418, t_419, pa_x, pc_x, fph0_416, fph0_417, fph0_418, \
                         fph0_419, fph1_416, fph1_417, fph1_418, \
                         fph1_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_416[k] = pa_x[k] * fph0_416[k]
                   - f_9 * pc_x[k] * fph1_416[k];

        t_417[k] = pa_x[k] * fph0_417[k]
                   - f_9 * pc_x[k] * fph1_417[k];

        t_418[k] = pa_x[k] * fph0_418[k]
                   - f_9 * pc_x[k] * fph1_418[k];

        t_419[k] = pa_x[k] * fph0_419[k]
                   - f_9 * pc_x[k] * fph1_419[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pb_z, pc_y, pc_z, fpg_165, gsh0_126, \
                         gsh0_129, gsg_90, gsh1_126, gsh1_129, \
                         gpg_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = pb_z[k] * gsh0_126[k]
                   - f_9 * pc_z[k] * gsh1_126[k];

        t_421[k] = f_11 * fpg_165[k]
                   + f_4 * pc_y[k] * gpg_300[k];

        t_422[k] = f_1 * gsg_90[k]
                   + f_4 * pc_z[k] * gpg_300[k];

        t_423[k] = pb_z[k] * gsh0_129[k]
                   - f_9 * pc_z[k] * gsh1_129[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, pa_x, pb_z, pc_x, pc_z, fph0_425, fpg_305, \
                         fph1_425, gsh0_132, gsg_91, gsh1_132, \
                         gpg_301 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = f_1 * gsg_91[k]
                   + f_4 * pc_z[k] * gpg_301[k];

        t_425[k] = pa_x[k] * fph0_425[k]
                   + f_11 * fpg_305[k]
                   - f_9 * pc_x[k] * fph1_425[k];

        t_426[k] = pb_z[k] * gsh0_132[k]
                   - f_9 * pc_z[k] * gsh1_132[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, pa_x, pc_x, pc_y, pc_z, fph0_429, fpg_170, \
                         fpg_309, fph1_429, gsg_93, gpg_303, gpg_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_1 * gsg_93[k]
                   + f_4 * pc_z[k] * gpg_303[k];

        t_428[k] = f_11 * fpg_170[k]
                   + f_4 * pc_y[k] * gpg_305[k];

        t_429[k] = pa_x[k] * fph0_429[k]
                   + f_10 * fpg_309[k]
                   - f_9 * pc_x[k] * fph1_429[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, pc_x, pc_z, fpg_310, fpg_312, fpg_313, \
                         gsg_96, gpg_306, gpg_310, gpg_312, gpg_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = f_1 * fpg_310[k]
                   + f_4 * pc_x[k] * gpg_310[k];

        t_431[k] = f_1 * gsg_96[k]
                   + f_4 * pc_z[k] * gpg_306[k];

        t_432[k] = f_1 * fpg_312[k]
                   + f_4 * pc_x[k] * gpg_312[k];

        t_433[k] = f_1 * fpg_313[k]
                   + f_4 * pc_x[k] * gpg_313[k];
    }

#pragma omp simd aligned(t_434, t_435, t_436, t_437, pa_x, pc_x, pc_z, fph0_435, fph0_437, \
                         fpg_314, fph1_435, fph1_437, gsg_100, gpg_310, \
                         gpg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_434[k] = f_1 * fpg_314[k]
                   + f_4 * pc_x[k] * gpg_314[k];

        t_435[k] = pa_x[k] * fph0_435[k]
                   - f_9 * pc_x[k] * fph1_435[k];

        t_436[k] = f_1 * gsg_100[k]
                   + f_4 * pc_z[k] * gpg_310[k];

        t_437[k] = pa_x[k] * fph0_437[k]
                   - f_9 * pc_x[k] * fph1_437[k];
    }

#pragma omp simd aligned(t_438, t_439, t_440, pa_x, pc_x, pc_y, fph0_438, fph0_440, fpg_179, \
                         fph1_438, fph1_440, gpg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_438[k] = pa_x[k] * fph0_438[k]
                   - f_9 * pc_x[k] * fph1_438[k];

        t_439[k] = f_11 * fpg_179[k]
                   + f_4 * pc_y[k] * gpg_314[k];

        t_440[k] = pa_x[k] * fph0_440[k]
                   - f_9 * pc_x[k] * fph1_440[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, pa_z, pc_y, pc_z, fph0_189, fph0_192, \
                         fpg_135, fpg_180, fph1_189, fph1_192, \
                         gpg_315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = pa_z[k] * fph0_189[k]
                   - f_9 * pc_z[k] * fph1_189[k];

        t_442[k] = f_10 * fpg_180[k]
                   + f_4 * pc_y[k] * gpg_315[k];

        t_443[k] = f_1 * fpg_135[k]
                   + f_4 * pc_z[k] * gpg_315[k];

        t_444[k] = pa_z[k] * fph0_192[k]
                   - f_9 * pc_z[k] * fph1_192[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, pa_z, pc_x, pc_y, pc_z, fph0_195, fpg_182, \
                         fpg_320, fph1_195, gsg_110, gpf0_215, gpf1_215, gpg_317, \
                         gpg_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_10 * fpg_182[k]
                   + f_4 * pc_y[k] * gpg_317[k];

        t_446[k] = f_1 * fpg_320[k]
                   + f_1 * gsg_110[k]
                   + f_7 * gpf0_215[k]
                   - f_8 * gpf1_215[k]
                   + f_4 * pc_x[k] * gpg_320[k];

        t_447[k] = pa_z[k] * fph0_195[k]
                   - f_9 * pc_z[k] * fph1_195[k];
    }

#pragma omp simd aligned(t_448, t_449, t_450, pc_x, pc_y, pc_z, fpg_138, fpg_185, fpg_324, \
                         gsg_114, gpf0_219, gpf1_219, gpg_318, gpg_320, \
                         gpg_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_448[k] = f_1 * fpg_138[k]
                   + f_4 * pc_z[k] * gpg_318[k];

        t_449[k] = f_10 * fpg_185[k]
                   + f_4 * pc_y[k] * gpg_320[k];

        t_450[k] = f_1 * fpg_324[k]
                   + f_1 * gsg_114[k]
                   + f_5 * gpf0_219[k]
                   - f_6 * gpf1_219[k]
                   + f_4 * pc_x[k] * gpg_324[k];
    }

#pragma omp simd aligned(t_451, t_452, t_453, pa_z, pc_x, pc_z, fph0_199, fpg_326, fpg_327, \
                         fph1_199, gsg_116, gsg_117, gpg_326, gpg_327 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_451[k] = pa_z[k] * fph0_199[k]
                   - f_9 * pc_z[k] * fph1_199[k];

        t_452[k] = f_1 * fpg_326[k]
                   + f_1 * gsg_116[k]
                   + f_4 * pc_x[k] * gpg_326[k];

        t_453[k] = f_1 * fpg_327[k]
                   + f_1 * gsg_117[k]
                   + f_4 * pc_x[k] * gpg_327[k];
    }

#pragma omp simd aligned(t_454, t_455, t_456, pc_x, pc_y, fpg_190, fpg_328, fpg_329, gsg_118, \
                         gsg_119, gpf0_216, gpf1_216, gpg_325, gpg_328, \
                         gpg_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_454[k] = f_1 * fpg_328[k]
                   + f_1 * gsg_118[k]
                   + f_4 * pc_x[k] * gpg_328[k];

        t_455[k] = f_1 * fpg_329[k]
                   + f_1 * gsg_119[k]
                   + f_4 * pc_x[k] * gpg_329[k];

        t_456[k] = f_10 * fpg_190[k]
                   + f_2 * gpf0_216[k]
                   - f_3 * gpf1_216[k]
                   + f_4 * pc_y[k] * gpg_325[k];
    }

#pragma omp simd aligned(t_457, t_458, t_459, pc_y, pc_z, fpg_145, fpg_192, fpg_193, gpf0_218, \
                         gpf0_219, gpf1_218, gpf1_219, gpg_325, gpg_327, \
                         gpg_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_457[k] = f_1 * fpg_145[k]
                   + f_4 * pc_z[k] * gpg_325[k];

        t_458[k] = f_10 * fpg_192[k]
                   + f_7 * gpf0_218[k]
                   - f_8 * gpf1_218[k]
                   + f_4 * pc_y[k] * gpg_327[k];

        t_459[k] = f_10 * fpg_193[k]
                   + f_5 * gpf0_219[k]
                   - f_6 * gpf1_219[k]
                   + f_4 * pc_y[k] * gpg_328[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, t_463, pa_z, pc_y, pc_z, fph0_210, fph0_211, \
                         fpg_149, fpg_194, fph1_210, fph1_211, gpf0_219, gpf1_219, \
                         gpg_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = f_10 * fpg_194[k]
                   + f_4 * pc_y[k] * gpg_329[k];

        t_461[k] = f_1 * fpg_149[k]
                   + f_2 * gpf0_219[k]
                   - f_3 * gpf1_219[k]
                   + f_4 * pc_z[k] * gpg_329[k];

        t_462[k] = pa_z[k] * fph0_210[k]
                   - f_9 * pc_z[k] * fph1_210[k];

        t_463[k] = pa_z[k] * fph0_211[k]
                   - f_9 * pc_z[k] * fph1_211[k];
    }

#pragma omp simd aligned(t_464, t_465, t_466, pa_z, pc_y, pc_z, fph0_213, fpg_150, fpg_197, \
                         fph1_213, gsg_107, gpg_330, gpg_332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_464[k] = f_1 * fpg_150[k]
                   + f_4 * pc_z[k] * gpg_330[k];

        t_465[k] = pa_z[k] * fph0_213[k]
                   - f_9 * pc_z[k] * fph1_213[k];

        t_466[k] = f_10 * fpg_197[k]
                   + f_1 * gsg_107[k]
                   + f_4 * pc_y[k] * gpg_332[k];
    }

#pragma omp simd aligned(t_467, t_468, t_469, pa_x, pa_z, pc_x, pc_z, fph0_216, fph0_467, \
                         fpg_153, fpg_335, fph1_216, fph1_467, \
                         gpg_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_467[k] = pa_x[k] * fph0_467[k]
                   + f_11 * fpg_335[k]
                   - f_9 * pc_x[k] * fph1_467[k];

        t_468[k] = pa_z[k] * fph0_216[k]
                   - f_9 * pc_z[k] * fph1_216[k];

        t_469[k] = f_1 * fpg_153[k]
                   + f_4 * pc_z[k] * gpg_333[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, pa_x, pc_x, pc_y, fph0_471, fpg_200, fpg_339, \
                         fpg_340, fph1_471, gsg_110, gpg_335, gpg_340 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_10 * fpg_200[k]
                   + f_1 * gsg_110[k]
                   + f_4 * pc_y[k] * gpg_335[k];

        t_471[k] = pa_x[k] * fph0_471[k]
                   + f_10 * fpg_339[k]
                   - f_9 * pc_x[k] * fph1_471[k];

        t_472[k] = f_1 * fpg_340[k]
                   + f_4 * pc_x[k] * gpg_340[k];
    }

#pragma omp simd aligned(t_473, t_474, t_475, t_476, pc_x, fpg_341, fpg_342, fpg_343, fpg_344, \
                         gpg_341, gpg_342, gpg_343, gpg_344 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_473[k] = f_1 * fpg_341[k]
                   + f_4 * pc_x[k] * gpg_341[k];

        t_474[k] = f_1 * fpg_342[k]
                   + f_4 * pc_x[k] * gpg_342[k];

        t_475[k] = f_1 * fpg_343[k]
                   + f_4 * pc_x[k] * gpg_343[k];

        t_476[k] = f_1 * fpg_344[k]
                   + f_4 * pc_x[k] * gpg_344[k];
    }
}

static auto
compute_prim_gph_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t fph0, const size_t fpg,
                                                          const size_t fph1, const size_t gsh0,
                                                          const size_t gsg, const size_t gsh1,
                                                          const size_t gpf0, const size_t gpf1,
                                                          const size_t gpg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 0.5 / q;
    const auto f_2 = 2.0 / gamma;
    const auto f_3 = 2.0 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = gamma / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_15 = 1.5 / gamma;
    const auto f_16 = 1.5 * p / (gamma * q);

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
    auto *t_592 = buffer.data(target + 592);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *fph0_315 = buffer.data(fph0 + 315);
    const auto *fph0_320 = buffer.data(fph0 + 320);
    const auto *fph0_324 = buffer.data(fph0 + 324);
    const auto *fph0_329 = buffer.data(fph0 + 329);
    const auto *fph0_357 = buffer.data(fph0 + 357);
    const auto *fph0_359 = buffer.data(fph0 + 359);
    const auto *fph0_362 = buffer.data(fph0 + 362);
    const auto *fph0_366 = buffer.data(fph0 + 366);
    const auto *fph0_477 = buffer.data(fph0 + 477);
    const auto *fph0_479 = buffer.data(fph0 + 479);
    const auto *fph0_480 = buffer.data(fph0 + 480);
    const auto *fph0_481 = buffer.data(fph0 + 481);
    const auto *fph0_482 = buffer.data(fph0 + 482);
    const auto *fph0_488 = buffer.data(fph0 + 488);
    const auto *fph0_492 = buffer.data(fph0 + 492);
    const auto *fph0_498 = buffer.data(fph0 + 498);
    const auto *fph0_499 = buffer.data(fph0 + 499);
    const auto *fph0_500 = buffer.data(fph0 + 500);
    const auto *fph0_501 = buffer.data(fph0 + 501);
    const auto *fph0_503 = buffer.data(fph0 + 503);
    const auto *fph0_528 = buffer.data(fph0 + 528);
    const auto *fph0_531 = buffer.data(fph0 + 531);
    const auto *fph0_540 = buffer.data(fph0 + 540);
    const auto *fph0_542 = buffer.data(fph0 + 542);
    const auto *fph0_543 = buffer.data(fph0 + 543);
    const auto *fph0_544 = buffer.data(fph0 + 544);
    const auto *fph0_545 = buffer.data(fph0 + 545);
    const auto *fph0_549 = buffer.data(fph0 + 549);
    const auto *fph0_552 = buffer.data(fph0 + 552);
    const auto *fph0_561 = buffer.data(fph0 + 561);
    const auto *fph0_562 = buffer.data(fph0 + 562);
    const auto *fph0_563 = buffer.data(fph0 + 563);
    const auto *fph0_564 = buffer.data(fph0 + 564);
    const auto *fph0_566 = buffer.data(fph0 + 566);
    const auto *fph0_591 = buffer.data(fph0 + 591);

    const auto *fpg_160 = buffer.data(fpg + 160);
    const auto *fpg_165 = buffer.data(fpg + 165);
    const auto *fpg_168 = buffer.data(fpg + 168);
    const auto *fpg_180 = buffer.data(fpg + 180);
    const auto *fpg_183 = buffer.data(fpg + 183);
    const auto *fpg_190 = buffer.data(fpg + 190);
    const auto *fpg_194 = buffer.data(fpg + 194);
    const auto *fpg_195 = buffer.data(fpg + 195);
    const auto *fpg_198 = buffer.data(fpg + 198);
    const auto *fpg_205 = buffer.data(fpg + 205);
    const auto *fpg_210 = buffer.data(fpg + 210);
    const auto *fpg_212 = buffer.data(fpg + 212);
    const auto *fpg_213 = buffer.data(fpg + 213);
    const auto *fpg_215 = buffer.data(fpg + 215);
    const auto *fpg_224 = buffer.data(fpg + 224);
    const auto *fpg_225 = buffer.data(fpg + 225);
    const auto *fpg_227 = buffer.data(fpg + 227);
    const auto *fpg_230 = buffer.data(fpg + 230);
    const auto *fpg_235 = buffer.data(fpg + 235);
    const auto *fpg_237 = buffer.data(fpg + 237);
    const auto *fpg_238 = buffer.data(fpg + 238);
    const auto *fpg_239 = buffer.data(fpg + 239);
    const auto *fpg_240 = buffer.data(fpg + 240);
    const auto *fpg_242 = buffer.data(fpg + 242);
    const auto *fpg_245 = buffer.data(fpg + 245);
    const auto *fpg_255 = buffer.data(fpg + 255);
    const auto *fpg_257 = buffer.data(fpg + 257);
    const auto *fpg_260 = buffer.data(fpg + 260);
    const auto *fpg_269 = buffer.data(fpg + 269);
    const auto *fpg_345 = buffer.data(fpg + 345);
    const auto *fpg_348 = buffer.data(fpg + 348);
    const auto *fpg_350 = buffer.data(fpg + 350);
    const auto *fpg_351 = buffer.data(fpg + 351);
    const auto *fpg_354 = buffer.data(fpg + 354);
    const auto *fpg_355 = buffer.data(fpg + 355);
    const auto *fpg_356 = buffer.data(fpg + 356);
    const auto *fpg_357 = buffer.data(fpg + 357);
    const auto *fpg_358 = buffer.data(fpg + 358);
    const auto *fpg_359 = buffer.data(fpg + 359);
    const auto *fpg_363 = buffer.data(fpg + 363);
    const auto *fpg_366 = buffer.data(fpg + 366);
    const auto *fpg_370 = buffer.data(fpg + 370);
    const auto *fpg_371 = buffer.data(fpg + 371);
    const auto *fpg_372 = buffer.data(fpg + 372);
    const auto *fpg_373 = buffer.data(fpg + 373);
    const auto *fpg_375 = buffer.data(fpg + 375);
    const auto *fpg_378 = buffer.data(fpg + 378);
    const auto *fpg_380 = buffer.data(fpg + 380);
    const auto *fpg_381 = buffer.data(fpg + 381);
    const auto *fpg_384 = buffer.data(fpg + 384);
    const auto *fpg_385 = buffer.data(fpg + 385);
    const auto *fpg_386 = buffer.data(fpg + 386);
    const auto *fpg_387 = buffer.data(fpg + 387);
    const auto *fpg_388 = buffer.data(fpg + 388);
    const auto *fpg_389 = buffer.data(fpg + 389);
    const auto *fpg_393 = buffer.data(fpg + 393);
    const auto *fpg_396 = buffer.data(fpg + 396);
    const auto *fpg_400 = buffer.data(fpg + 400);
    const auto *fpg_401 = buffer.data(fpg + 401);
    const auto *fpg_402 = buffer.data(fpg + 402);
    const auto *fpg_403 = buffer.data(fpg + 403);
    const auto *fpg_404 = buffer.data(fpg + 404);
    const auto *fpg_405 = buffer.data(fpg + 405);
    const auto *fpg_410 = buffer.data(fpg + 410);
    const auto *fpg_414 = buffer.data(fpg + 414);
    const auto *fpg_415 = buffer.data(fpg + 415);
    const auto *fpg_416 = buffer.data(fpg + 416);
    const auto *fpg_417 = buffer.data(fpg + 417);
    const auto *fpg_419 = buffer.data(fpg + 419);
    const auto *fpg_423 = buffer.data(fpg + 423);

    const auto *fph1_315 = buffer.data(fph1 + 315);
    const auto *fph1_320 = buffer.data(fph1 + 320);
    const auto *fph1_324 = buffer.data(fph1 + 324);
    const auto *fph1_329 = buffer.data(fph1 + 329);
    const auto *fph1_357 = buffer.data(fph1 + 357);
    const auto *fph1_359 = buffer.data(fph1 + 359);
    const auto *fph1_362 = buffer.data(fph1 + 362);
    const auto *fph1_366 = buffer.data(fph1 + 366);
    const auto *fph1_477 = buffer.data(fph1 + 477);
    const auto *fph1_479 = buffer.data(fph1 + 479);
    const auto *fph1_480 = buffer.data(fph1 + 480);
    const auto *fph1_481 = buffer.data(fph1 + 481);
    const auto *fph1_482 = buffer.data(fph1 + 482);
    const auto *fph1_488 = buffer.data(fph1 + 488);
    const auto *fph1_492 = buffer.data(fph1 + 492);
    const auto *fph1_498 = buffer.data(fph1 + 498);
    const auto *fph1_499 = buffer.data(fph1 + 499);
    const auto *fph1_500 = buffer.data(fph1 + 500);
    const auto *fph1_501 = buffer.data(fph1 + 501);
    const auto *fph1_503 = buffer.data(fph1 + 503);
    const auto *fph1_528 = buffer.data(fph1 + 528);
    const auto *fph1_531 = buffer.data(fph1 + 531);
    const auto *fph1_540 = buffer.data(fph1 + 540);
    const auto *fph1_542 = buffer.data(fph1 + 542);
    const auto *fph1_543 = buffer.data(fph1 + 543);
    const auto *fph1_544 = buffer.data(fph1 + 544);
    const auto *fph1_545 = buffer.data(fph1 + 545);
    const auto *fph1_549 = buffer.data(fph1 + 549);
    const auto *fph1_552 = buffer.data(fph1 + 552);
    const auto *fph1_561 = buffer.data(fph1 + 561);
    const auto *fph1_562 = buffer.data(fph1 + 562);
    const auto *fph1_563 = buffer.data(fph1 + 563);
    const auto *fph1_564 = buffer.data(fph1 + 564);
    const auto *fph1_566 = buffer.data(fph1 + 566);
    const auto *fph1_591 = buffer.data(fph1 + 591);

    const auto *gsh0_189 = buffer.data(gsh0 + 189);

    const auto *gsg_105 = buffer.data(gsg + 105);
    const auto *gsg_108 = buffer.data(gsg + 108);
    const auto *gsg_120 = buffer.data(gsg + 120);
    const auto *gsg_122 = buffer.data(gsg + 122);
    const auto *gsg_123 = buffer.data(gsg + 123);
    const auto *gsg_125 = buffer.data(gsg + 125);
    const auto *gsg_126 = buffer.data(gsg + 126);
    const auto *gsg_130 = buffer.data(gsg + 130);
    const auto *gsg_131 = buffer.data(gsg + 131);
    const auto *gsg_132 = buffer.data(gsg + 132);
    const auto *gsg_133 = buffer.data(gsg + 133);
    const auto *gsg_135 = buffer.data(gsg + 135);
    const auto *gsg_137 = buffer.data(gsg + 137);
    const auto *gsg_140 = buffer.data(gsg + 140);
    const auto *gsg_144 = buffer.data(gsg + 144);
    const auto *gsg_145 = buffer.data(gsg + 145);
    const auto *gsg_146 = buffer.data(gsg + 146);
    const auto *gsg_147 = buffer.data(gsg + 147);
    const auto *gsg_149 = buffer.data(gsg + 149);

    const auto *gsh1_189 = buffer.data(gsh1 + 189);

    const auto *gpf0_230 = buffer.data(gpf0 + 230);
    const auto *gpf0_233 = buffer.data(gpf0 + 233);
    const auto *gpf0_236 = buffer.data(gpf0 + 236);
    const auto *gpf0_243 = buffer.data(gpf0 + 243);
    const auto *gpf0_246 = buffer.data(gpf0 + 246);
    const auto *gpf0_248 = buffer.data(gpf0 + 248);
    const auto *gpf0_249 = buffer.data(gpf0 + 249);
    const auto *gpf0_250 = buffer.data(gpf0 + 250);
    const auto *gpf0_255 = buffer.data(gpf0 + 255);
    const auto *gpf0_259 = buffer.data(gpf0 + 259);
    const auto *gpf0_270 = buffer.data(gpf0 + 270);
    const auto *gpf0_271 = buffer.data(gpf0 + 271);
    const auto *gpf0_272 = buffer.data(gpf0 + 272);
    const auto *gpf0_275 = buffer.data(gpf0 + 275);
    const auto *gpf0_276 = buffer.data(gpf0 + 276);
    const auto *gpf0_277 = buffer.data(gpf0 + 277);
    const auto *gpf0_278 = buffer.data(gpf0 + 278);
    const auto *gpf0_279 = buffer.data(gpf0 + 279);

    const auto *gpf1_230 = buffer.data(gpf1 + 230);
    const auto *gpf1_233 = buffer.data(gpf1 + 233);
    const auto *gpf1_236 = buffer.data(gpf1 + 236);
    const auto *gpf1_243 = buffer.data(gpf1 + 243);
    const auto *gpf1_246 = buffer.data(gpf1 + 246);
    const auto *gpf1_248 = buffer.data(gpf1 + 248);
    const auto *gpf1_249 = buffer.data(gpf1 + 249);
    const auto *gpf1_250 = buffer.data(gpf1 + 250);
    const auto *gpf1_255 = buffer.data(gpf1 + 255);
    const auto *gpf1_259 = buffer.data(gpf1 + 259);
    const auto *gpf1_270 = buffer.data(gpf1 + 270);
    const auto *gpf1_271 = buffer.data(gpf1 + 271);
    const auto *gpf1_272 = buffer.data(gpf1 + 272);
    const auto *gpf1_275 = buffer.data(gpf1 + 275);
    const auto *gpf1_276 = buffer.data(gpf1 + 276);
    const auto *gpf1_277 = buffer.data(gpf1 + 277);
    const auto *gpf1_278 = buffer.data(gpf1 + 278);
    const auto *gpf1_279 = buffer.data(gpf1 + 279);

    const auto *gpg_340 = buffer.data(gpg + 340);
    const auto *gpg_345 = buffer.data(gpg + 345);
    const auto *gpg_347 = buffer.data(gpg + 347);
    const auto *gpg_348 = buffer.data(gpg + 348);
    const auto *gpg_350 = buffer.data(gpg + 350);
    const auto *gpg_351 = buffer.data(gpg + 351);
    const auto *gpg_355 = buffer.data(gpg + 355);
    const auto *gpg_356 = buffer.data(gpg + 356);
    const auto *gpg_357 = buffer.data(gpg + 357);
    const auto *gpg_358 = buffer.data(gpg + 358);
    const auto *gpg_359 = buffer.data(gpg + 359);
    const auto *gpg_360 = buffer.data(gpg + 360);
    const auto *gpg_362 = buffer.data(gpg + 362);
    const auto *gpg_363 = buffer.data(gpg + 363);
    const auto *gpg_365 = buffer.data(gpg + 365);
    const auto *gpg_366 = buffer.data(gpg + 366);
    const auto *gpg_370 = buffer.data(gpg + 370);
    const auto *gpg_371 = buffer.data(gpg + 371);
    const auto *gpg_372 = buffer.data(gpg + 372);
    const auto *gpg_373 = buffer.data(gpg + 373);
    const auto *gpg_374 = buffer.data(gpg + 374);
    const auto *gpg_375 = buffer.data(gpg + 375);
    const auto *gpg_377 = buffer.data(gpg + 377);
    const auto *gpg_378 = buffer.data(gpg + 378);
    const auto *gpg_380 = buffer.data(gpg + 380);
    const auto *gpg_384 = buffer.data(gpg + 384);
    const auto *gpg_385 = buffer.data(gpg + 385);
    const auto *gpg_386 = buffer.data(gpg + 386);
    const auto *gpg_387 = buffer.data(gpg + 387);
    const auto *gpg_388 = buffer.data(gpg + 388);
    const auto *gpg_389 = buffer.data(gpg + 389);
    const auto *gpg_390 = buffer.data(gpg + 390);
    const auto *gpg_392 = buffer.data(gpg + 392);
    const auto *gpg_393 = buffer.data(gpg + 393);
    const auto *gpg_395 = buffer.data(gpg + 395);
    const auto *gpg_400 = buffer.data(gpg + 400);
    const auto *gpg_401 = buffer.data(gpg + 401);
    const auto *gpg_402 = buffer.data(gpg + 402);
    const auto *gpg_403 = buffer.data(gpg + 403);
    const auto *gpg_404 = buffer.data(gpg + 404);
    const auto *gpg_405 = buffer.data(gpg + 405);
    const auto *gpg_406 = buffer.data(gpg + 406);
    const auto *gpg_407 = buffer.data(gpg + 407);
    const auto *gpg_408 = buffer.data(gpg + 408);
    const auto *gpg_409 = buffer.data(gpg + 409);
    const auto *gpg_410 = buffer.data(gpg + 410);
    const auto *gpg_414 = buffer.data(gpg + 414);
    const auto *gpg_415 = buffer.data(gpg + 415);
    const auto *gpg_416 = buffer.data(gpg + 416);
    const auto *gpg_417 = buffer.data(gpg + 417);
    const auto *gpg_418 = buffer.data(gpg + 418);
    const auto *gpg_419 = buffer.data(gpg + 419);
    const auto *gpg_420 = buffer.data(gpg + 420);
    const auto *gpg_422 = buffer.data(gpg + 422);

#pragma omp simd aligned(t_477, t_478, t_479, t_480, pa_x, pc_x, pc_z, fph0_477, fph0_479, \
                         fph0_480, fpg_160, fph1_477, fph1_479, fph1_480, \
                         gpg_340 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_477[k] = pa_x[k] * fph0_477[k]
                   - f_9 * pc_x[k] * fph1_477[k];

        t_478[k] = f_1 * fpg_160[k]
                   + f_4 * pc_z[k] * gpg_340[k];

        t_479[k] = pa_x[k] * fph0_479[k]
                   - f_9 * pc_x[k] * fph1_479[k];

        t_480[k] = pa_x[k] * fph0_480[k]
                   - f_9 * pc_x[k] * fph1_480[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, t_484, pa_x, pc_x, pc_y, fph0_481, fph0_482, \
                         fpg_210, fpg_345, fph1_481, fph1_482, gpf0_230, gpf1_230, \
                         gpg_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = pa_x[k] * fph0_481[k]
                   - f_9 * pc_x[k] * fph1_481[k];

        t_482[k] = pa_x[k] * fph0_482[k]
                   - f_9 * pc_x[k] * fph1_482[k];

        t_483[k] = f_1 * fpg_345[k]
                   + f_2 * gpf0_230[k]
                   - f_3 * gpf1_230[k]
                   + f_4 * pc_x[k] * gpg_345[k];

        t_484[k] = f_10 * fpg_210[k]
                   + f_4 * pc_y[k] * gpg_345[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, pc_x, pc_y, pc_z, fpg_165, fpg_212, fpg_348, \
                         gsg_105, gpf0_233, gpf1_233, gpg_345, gpg_347, \
                         gpg_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_1 * fpg_165[k]
                   + f_1 * gsg_105[k]
                   + f_4 * pc_z[k] * gpg_345[k];

        t_486[k] = f_1 * fpg_348[k]
                   + f_7 * gpf0_233[k]
                   - f_8 * gpf1_233[k]
                   + f_4 * pc_x[k] * gpg_348[k];

        t_487[k] = f_10 * fpg_212[k]
                   + f_4 * pc_y[k] * gpg_347[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, pa_x, pc_x, pc_z, fph0_488, fpg_168, fpg_350, \
                         fpg_351, fph1_488, gsg_108, gpf0_236, gpf1_236, gpg_348, \
                         gpg_351 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = pa_x[k] * fph0_488[k]
                   + f_11 * fpg_350[k]
                   - f_9 * pc_x[k] * fph1_488[k];

        t_489[k] = f_1 * fpg_351[k]
                   + f_5 * gpf0_236[k]
                   - f_6 * gpf1_236[k]
                   + f_4 * pc_x[k] * gpg_351[k];

        t_490[k] = f_1 * fpg_168[k]
                   + f_1 * gsg_108[k]
                   + f_4 * pc_z[k] * gpg_348[k];
    }

#pragma omp simd aligned(t_491, t_492, t_493, t_494, pa_x, pc_x, pc_y, fph0_492, fpg_215, \
                         fpg_354, fpg_355, fpg_356, fph1_492, gpg_350, gpg_355, \
                         gpg_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_491[k] = f_10 * fpg_215[k]
                   + f_4 * pc_y[k] * gpg_350[k];

        t_492[k] = pa_x[k] * fph0_492[k]
                   + f_10 * fpg_354[k]
                   - f_9 * pc_x[k] * fph1_492[k];

        t_493[k] = f_1 * fpg_355[k]
                   + f_4 * pc_x[k] * gpg_355[k];

        t_494[k] = f_1 * fpg_356[k]
                   + f_4 * pc_x[k] * gpg_356[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, t_498, pa_x, pc_x, fph0_498, fpg_357, fpg_358, \
                         fpg_359, fph1_498, gpg_357, gpg_358, gpg_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = f_1 * fpg_357[k]
                   + f_4 * pc_x[k] * gpg_357[k];

        t_496[k] = f_1 * fpg_358[k]
                   + f_4 * pc_x[k] * gpg_358[k];

        t_497[k] = f_1 * fpg_359[k]
                   + f_4 * pc_x[k] * gpg_359[k];

        t_498[k] = pa_x[k] * fph0_498[k]
                   - f_9 * pc_x[k] * fph1_498[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, t_502, pa_x, pc_x, pc_y, fph0_499, fph0_500, \
                         fph0_501, fpg_224, fph1_499, fph1_500, fph1_501, \
                         gpg_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = pa_x[k] * fph0_499[k]
                   - f_9 * pc_x[k] * fph1_499[k];

        t_500[k] = pa_x[k] * fph0_500[k]
                   - f_9 * pc_x[k] * fph1_500[k];

        t_501[k] = pa_x[k] * fph0_501[k]
                   - f_9 * pc_x[k] * fph1_501[k];

        t_502[k] = f_10 * fpg_224[k]
                   + f_4 * pc_y[k] * gpg_359[k];
    }

#pragma omp simd aligned(t_503, t_504, t_505, t_506, pa_x, pa_y, pc_x, pc_y, pc_z, fph0_315, \
                         fph0_503, fpg_180, fpg_225, fph1_315, fph1_503, \
                         gpg_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_503[k] = pa_x[k] * fph0_503[k]
                   - f_9 * pc_x[k] * fph1_503[k];

        t_504[k] = pa_y[k] * fph0_315[k]
                   - f_9 * pc_y[k] * fph1_315[k];

        t_505[k] = f_1 * fpg_225[k]
                   + f_4 * pc_y[k] * gpg_360[k];

        t_506[k] = f_10 * fpg_180[k]
                   + f_4 * pc_z[k] * gpg_360[k];
    }

#pragma omp simd aligned(t_507, t_508, t_509, pa_y, pc_x, pc_y, fph0_320, fpg_227, fpg_363, \
                         fph1_320, gsg_123, gpf0_243, gpf1_243, gpg_362, \
                         gpg_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_507[k] = f_1 * fpg_363[k]
                   + f_1 * gsg_123[k]
                   + f_7 * gpf0_243[k]
                   - f_8 * gpf1_243[k]
                   + f_4 * pc_x[k] * gpg_363[k];

        t_508[k] = f_1 * fpg_227[k]
                   + f_4 * pc_y[k] * gpg_362[k];

        t_509[k] = pa_y[k] * fph0_320[k]
                   - f_9 * pc_y[k] * fph1_320[k];
    }

#pragma omp simd aligned(t_510, t_511, t_512, pc_x, pc_y, pc_z, fpg_183, fpg_230, fpg_366, \
                         gsg_126, gpf0_246, gpf1_246, gpg_363, gpg_365, \
                         gpg_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_510[k] = f_1 * fpg_366[k]
                   + f_1 * gsg_126[k]
                   + f_5 * gpf0_246[k]
                   - f_6 * gpf1_246[k]
                   + f_4 * pc_x[k] * gpg_366[k];

        t_511[k] = f_10 * fpg_183[k]
                   + f_4 * pc_z[k] * gpg_363[k];

        t_512[k] = f_1 * fpg_230[k]
                   + f_4 * pc_y[k] * gpg_365[k];
    }

#pragma omp simd aligned(t_513, t_514, t_515, pa_y, pc_x, pc_y, fph0_324, fpg_370, fpg_371, \
                         fph1_324, gsg_130, gsg_131, gpg_370, gpg_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_513[k] = pa_y[k] * fph0_324[k]
                   - f_9 * pc_y[k] * fph1_324[k];

        t_514[k] = f_1 * fpg_370[k]
                   + f_1 * gsg_130[k]
                   + f_4 * pc_x[k] * gpg_370[k];

        t_515[k] = f_1 * fpg_371[k]
                   + f_1 * gsg_131[k]
                   + f_4 * pc_x[k] * gpg_371[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, pa_y, pc_x, pc_y, fph0_329, fpg_372, fpg_373, \
                         fph1_329, gsg_132, gsg_133, gpg_372, gpg_373 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = f_1 * fpg_372[k]
                   + f_1 * gsg_132[k]
                   + f_4 * pc_x[k] * gpg_372[k];

        t_517[k] = f_1 * fpg_373[k]
                   + f_1 * gsg_133[k]
                   + f_4 * pc_x[k] * gpg_373[k];

        t_518[k] = pa_y[k] * fph0_329[k]
                   - f_9 * pc_y[k] * fph1_329[k];
    }

#pragma omp simd aligned(t_519, t_520, t_521, pc_y, pc_z, fpg_190, fpg_235, fpg_237, gpf0_246, \
                         gpf0_248, gpf1_246, gpf1_248, gpg_370, \
                         gpg_372 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_519[k] = f_1 * fpg_235[k]
                   + f_2 * gpf0_246[k]
                   - f_3 * gpf1_246[k]
                   + f_4 * pc_y[k] * gpg_370[k];

        t_520[k] = f_10 * fpg_190[k]
                   + f_4 * pc_z[k] * gpg_370[k];

        t_521[k] = f_1 * fpg_237[k]
                   + f_7 * gpf0_248[k]
                   - f_8 * gpf1_248[k]
                   + f_4 * pc_y[k] * gpg_372[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, pc_y, pc_z, fpg_194, fpg_238, fpg_239, gpf0_249, \
                         gpf1_249, gpg_373, gpg_374 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_1 * fpg_238[k]
                   + f_5 * gpf0_249[k]
                   - f_6 * gpf1_249[k]
                   + f_4 * pc_y[k] * gpg_373[k];

        t_523[k] = f_1 * fpg_239[k]
                   + f_4 * pc_y[k] * gpg_374[k];

        t_524[k] = f_10 * fpg_194[k]
                   + f_2 * gpf0_249[k]
                   - f_3 * gpf1_249[k]
                   + f_4 * pc_z[k] * gpg_374[k];
    }

#pragma omp simd aligned(t_525, t_526, t_527, pc_x, pc_y, pc_z, fpg_195, fpg_240, fpg_375, \
                         gsg_120, gpf0_250, gpf1_250, gpg_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_525[k] = f_1 * fpg_375[k]
                   + f_2 * gpf0_250[k]
                   - f_3 * gpf1_250[k]
                   + f_4 * pc_x[k] * gpg_375[k];

        t_526[k] = f_1 * fpg_240[k]
                   + f_1 * gsg_120[k]
                   + f_4 * pc_y[k] * gpg_375[k];

        t_527[k] = f_10 * fpg_195[k]
                   + f_4 * pc_z[k] * gpg_375[k];
    }

#pragma omp simd aligned(t_528, t_529, t_530, pa_x, pc_x, pc_y, fph0_528, fpg_242, fpg_378, \
                         fpg_380, fph1_528, gsg_122, gpf0_255, gpf1_255, gpg_377, \
                         gpg_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_528[k] = pa_x[k] * fph0_528[k]
                   + f_11 * fpg_378[k]
                   - f_9 * pc_x[k] * fph1_528[k];

        t_529[k] = f_1 * fpg_242[k]
                   + f_1 * gsg_122[k]
                   + f_4 * pc_y[k] * gpg_377[k];

        t_530[k] = f_1 * fpg_380[k]
                   + f_7 * gpf0_255[k]
                   - f_8 * gpf1_255[k]
                   + f_4 * pc_x[k] * gpg_380[k];
    }

#pragma omp simd aligned(t_531, t_532, t_533, pa_x, pc_x, pc_y, pc_z, fph0_531, fpg_198, \
                         fpg_245, fpg_381, fph1_531, gsg_125, gpg_378, \
                         gpg_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_531[k] = pa_x[k] * fph0_531[k]
                   + f_10 * fpg_381[k]
                   - f_9 * pc_x[k] * fph1_531[k];

        t_532[k] = f_10 * fpg_198[k]
                   + f_4 * pc_z[k] * gpg_378[k];

        t_533[k] = f_1 * fpg_245[k]
                   + f_1 * gsg_125[k]
                   + f_4 * pc_y[k] * gpg_380[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, t_537, pc_x, fpg_384, fpg_385, fpg_386, fpg_387, \
                         gpf0_259, gpf1_259, gpg_384, gpg_385, gpg_386, \
                         gpg_387 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = f_1 * fpg_384[k]
                   + f_5 * gpf0_259[k]
                   - f_6 * gpf1_259[k]
                   + f_4 * pc_x[k] * gpg_384[k];

        t_535[k] = f_1 * fpg_385[k]
                   + f_4 * pc_x[k] * gpg_385[k];

        t_536[k] = f_1 * fpg_386[k]
                   + f_4 * pc_x[k] * gpg_386[k];

        t_537[k] = f_1 * fpg_387[k]
                   + f_4 * pc_x[k] * gpg_387[k];
    }

#pragma omp simd aligned(t_538, t_539, t_540, t_541, pa_x, pc_x, pc_z, fph0_540, fpg_205, \
                         fpg_388, fpg_389, fph1_540, gpg_385, gpg_388, \
                         gpg_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_538[k] = f_1 * fpg_388[k]
                   + f_4 * pc_x[k] * gpg_388[k];

        t_539[k] = f_1 * fpg_389[k]
                   + f_4 * pc_x[k] * gpg_389[k];

        t_540[k] = pa_x[k] * fph0_540[k]
                   - f_9 * pc_x[k] * fph1_540[k];

        t_541[k] = f_10 * fpg_205[k]
                   + f_4 * pc_z[k] * gpg_385[k];
    }

#pragma omp simd aligned(t_542, t_543, t_544, t_545, pa_x, pc_x, fph0_542, fph0_543, fph0_544, \
                         fph0_545, fph1_542, fph1_543, fph1_544, \
                         fph1_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = pa_x[k] * fph0_542[k]
                   - f_9 * pc_x[k] * fph1_542[k];

        t_543[k] = pa_x[k] * fph0_543[k]
                   - f_9 * pc_x[k] * fph1_543[k];

        t_544[k] = pa_x[k] * fph0_544[k]
                   - f_9 * pc_x[k] * fph1_544[k];

        t_545[k] = pa_x[k] * fph0_545[k]
                   - f_9 * pc_x[k] * fph1_545[k];
    }

#pragma omp simd aligned(t_546, t_547, t_548, pa_y, pc_y, fph0_357, fph0_359, fpg_255, \
                         fph1_357, fph1_359, gpg_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_546[k] = pa_y[k] * fph0_357[k]
                   - f_9 * pc_y[k] * fph1_357[k];

        t_547[k] = f_1 * fpg_255[k]
                   + f_4 * pc_y[k] * gpg_390[k];

        t_548[k] = pa_y[k] * fph0_359[k]
                   - f_9 * pc_y[k] * fph1_359[k];
    }

#pragma omp simd aligned(t_549, t_550, t_551, pa_x, pa_y, pc_x, pc_y, fph0_362, fph0_549, \
                         fpg_257, fpg_393, fph1_362, fph1_549, \
                         gpg_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_549[k] = pa_x[k] * fph0_549[k]
                   + f_11 * fpg_393[k]
                   - f_9 * pc_x[k] * fph1_549[k];

        t_550[k] = f_1 * fpg_257[k]
                   + f_4 * pc_y[k] * gpg_392[k];

        t_551[k] = pa_y[k] * fph0_362[k]
                   - f_9 * pc_y[k] * fph1_362[k];
    }

#pragma omp simd aligned(t_552, t_553, t_554, pa_x, pc_x, pc_y, pc_z, fph0_552, fpg_213, \
                         fpg_260, fpg_396, fph1_552, gsg_123, gpg_393, \
                         gpg_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_552[k] = pa_x[k] * fph0_552[k]
                   + f_10 * fpg_396[k]
                   - f_9 * pc_x[k] * fph1_552[k];

        t_553[k] = f_10 * fpg_213[k]
                   + f_1 * gsg_123[k]
                   + f_4 * pc_z[k] * gpg_393[k];

        t_554[k] = f_1 * fpg_260[k]
                   + f_4 * pc_y[k] * gpg_395[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, t_558, pa_y, pc_x, pc_y, fph0_366, fpg_400, \
                         fpg_401, fpg_402, fph1_366, gpg_400, gpg_401, \
                         gpg_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = pa_y[k] * fph0_366[k]
                   - f_9 * pc_y[k] * fph1_366[k];

        t_556[k] = f_1 * fpg_400[k]
                   + f_4 * pc_x[k] * gpg_400[k];

        t_557[k] = f_1 * fpg_401[k]
                   + f_4 * pc_x[k] * gpg_401[k];

        t_558[k] = f_1 * fpg_402[k]
                   + f_4 * pc_x[k] * gpg_402[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, t_562, pa_x, pc_x, fph0_561, fph0_562, fpg_403, \
                         fpg_404, fph1_561, fph1_562, gpg_403, \
                         gpg_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = f_1 * fpg_403[k]
                   + f_4 * pc_x[k] * gpg_403[k];

        t_560[k] = f_1 * fpg_404[k]
                   + f_4 * pc_x[k] * gpg_404[k];

        t_561[k] = pa_x[k] * fph0_561[k]
                   - f_9 * pc_x[k] * fph1_561[k];

        t_562[k] = pa_x[k] * fph0_562[k]
                   - f_9 * pc_x[k] * fph1_562[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, t_566, pa_x, pc_x, pc_y, fph0_563, fph0_564, \
                         fph0_566, fpg_269, fph1_563, fph1_564, fph1_566, \
                         gpg_404 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = pa_x[k] * fph0_563[k]
                   - f_9 * pc_x[k] * fph1_563[k];

        t_564[k] = pa_x[k] * fph0_564[k]
                   - f_9 * pc_x[k] * fph1_564[k];

        t_565[k] = f_1 * fpg_269[k]
                   + f_4 * pc_y[k] * gpg_404[k];

        t_566[k] = pa_x[k] * fph0_566[k]
                   - f_9 * pc_x[k] * fph1_566[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, t_570, t_571, pc_x, pc_y, pc_z, fpg_225, \
                         fpg_405, gsg_135, gpf0_270, gpf1_270, gpg_405, gpg_406, \
                         gpg_407 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = f_1 * fpg_405[k]
                   + f_1 * gsg_135[k]
                   + f_2 * gpf0_270[k]
                   - f_3 * gpf1_270[k]
                   + f_4 * pc_x[k] * gpg_405[k];

        t_568[k] = f_4 * pc_y[k] * gpg_405[k];

        t_569[k] = f_11 * fpg_225[k]
                   + f_4 * pc_z[k] * gpg_405[k];

        t_570[k] = f_5 * gpf0_270[k]
                   - f_6 * gpf1_270[k]
                   + f_4 * pc_y[k] * gpg_406[k];

        t_571[k] = f_4 * pc_y[k] * gpg_407[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, pc_x, pc_y, fpg_410, gsg_140, gpf0_271, \
                         gpf0_272, gpf0_275, gpf1_271, gpf1_272, gpf1_275, gpg_408, gpg_409, \
                         gpg_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = f_1 * fpg_410[k]
                   + f_1 * gsg_140[k]
                   + f_7 * gpf0_275[k]
                   - f_8 * gpf1_275[k]
                   + f_4 * pc_x[k] * gpg_410[k];

        t_573[k] = f_7 * gpf0_271[k]
                   - f_8 * gpf1_271[k]
                   + f_4 * pc_y[k] * gpg_408[k];

        t_574[k] = f_5 * gpf0_272[k]
                   - f_6 * gpf1_272[k]
                   + f_4 * pc_y[k] * gpg_409[k];
    }

#pragma omp simd aligned(t_575, t_576, t_577, pc_x, pc_y, fpg_414, fpg_415, gsg_144, gsg_145, \
                         gpf0_279, gpf1_279, gpg_410, gpg_414, \
                         gpg_415 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = f_4 * pc_y[k] * gpg_410[k];

        t_576[k] = f_1 * fpg_414[k]
                   + f_1 * gsg_144[k]
                   + f_5 * gpf0_279[k]
                   - f_6 * gpf1_279[k]
                   + f_4 * pc_x[k] * gpg_414[k];

        t_577[k] = f_1 * fpg_415[k]
                   + f_1 * gsg_145[k]
                   + f_4 * pc_x[k] * gpg_415[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, t_581, pc_x, pc_y, fpg_416, fpg_417, fpg_419, \
                         gsg_146, gsg_147, gsg_149, gpg_414, gpg_416, gpg_417, \
                         gpg_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_1 * fpg_416[k]
                   + f_1 * gsg_146[k]
                   + f_4 * pc_x[k] * gpg_416[k];

        t_579[k] = f_1 * fpg_417[k]
                   + f_1 * gsg_147[k]
                   + f_4 * pc_x[k] * gpg_417[k];

        t_580[k] = f_4 * pc_y[k] * gpg_414[k];

        t_581[k] = f_1 * fpg_419[k]
                   + f_1 * gsg_149[k]
                   + f_4 * pc_x[k] * gpg_419[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, pc_y, gpf0_276, gpf0_277, gpf0_278, gpf1_276, \
                         gpf1_277, gpf1_278, gpg_415, gpg_416, \
                         gpg_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = f_2 * gpf0_276[k]
                   - f_3 * gpf1_276[k]
                   + f_4 * pc_y[k] * gpg_415[k];

        t_583[k] = f_15 * gpf0_277[k]
                   - f_16 * gpf1_277[k]
                   + f_4 * pc_y[k] * gpg_416[k];

        t_584[k] = f_7 * gpf0_278[k]
                   - f_8 * gpf1_278[k]
                   + f_4 * pc_y[k] * gpg_417[k];
    }

#pragma omp simd aligned(t_585, t_586, t_587, t_588, pb_y, pc_y, pc_z, fpg_239, gsh0_189, \
                         gsh1_189, gpf0_279, gpf1_279, gpg_418, \
                         gpg_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_585[k] = f_5 * gpf0_279[k]
                   - f_6 * gpf1_279[k]
                   + f_4 * pc_y[k] * gpg_418[k];

        t_586[k] = f_4 * pc_y[k] * gpg_419[k];

        t_587[k] = f_11 * fpg_239[k]
                   + f_2 * gpf0_279[k]
                   - f_3 * gpf1_279[k]
                   + f_4 * pc_z[k] * gpg_419[k];

        t_588[k] = pb_y[k] * gsh0_189[k]
                   - f_9 * pc_y[k] * gsh1_189[k];
    }

#pragma omp simd aligned(t_589, t_590, t_591, t_592, pa_x, pc_x, pc_y, pc_z, fph0_591, \
                         fpg_240, fpg_423, fph1_591, gsg_135, gsg_137, gpg_420, \
                         gpg_422 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_589[k] = f_1 * gsg_135[k]
                   + f_4 * pc_y[k] * gpg_420[k];

        t_590[k] = f_11 * fpg_240[k]
                   + f_4 * pc_z[k] * gpg_420[k];

        t_591[k] = pa_x[k] * fph0_591[k]
                   + f_11 * fpg_423[k]
                   - f_9 * pc_x[k] * fph1_591[k];

        t_592[k] = f_1 * gsg_137[k]
                   + f_4 * pc_y[k] * gpg_422[k];
    }
}

static auto
compute_prim_gph_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t fph0, const size_t fpg,
                                                          const size_t fph1, const size_t gsh0,
                                                          const size_t gsg, const size_t gsh1,
                                                          const size_t gpf0, const size_t gpf1,
                                                          const size_t gpg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 2.0 / gamma;
    const auto f_3 = 2.0 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = gamma / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 2.5 / q;
    const auto f_15 = 1.5 / gamma;
    const auto f_16 = 1.5 * p / (gamma * q);

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
    auto *t_716 = buffer.data(target + 716);
    auto *t_717 = buffer.data(target + 717);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *fph0_378 = buffer.data(fph0 + 378);
    const auto *fph0_379 = buffer.data(fph0 + 379);
    const auto *fph0_381 = buffer.data(fph0 + 381);
    const auto *fph0_382 = buffer.data(fph0 + 382);
    const auto *fph0_384 = buffer.data(fph0 + 384);
    const auto *fph0_385 = buffer.data(fph0 + 385);
    const auto *fph0_386 = buffer.data(fph0 + 386);
    const auto *fph0_393 = buffer.data(fph0 + 393);
    const auto *fph0_399 = buffer.data(fph0 + 399);
    const auto *fph0_400 = buffer.data(fph0 + 400);
    const auto *fph0_402 = buffer.data(fph0 + 402);
    const auto *fph0_594 = buffer.data(fph0 + 594);
    const auto *fph0_595 = buffer.data(fph0 + 595);
    const auto *fph0_603 = buffer.data(fph0 + 603);
    const auto *fph0_604 = buffer.data(fph0 + 604);
    const auto *fph0_605 = buffer.data(fph0 + 605);
    const auto *fph0_606 = buffer.data(fph0 + 606);
    const auto *fph0_608 = buffer.data(fph0 + 608);
    const auto *fph0_609 = buffer.data(fph0 + 609);
    const auto *fph0_614 = buffer.data(fph0 + 614);
    const auto *fph0_618 = buffer.data(fph0 + 618);
    const auto *fph0_624 = buffer.data(fph0 + 624);
    const auto *fph0_625 = buffer.data(fph0 + 625);
    const auto *fph0_626 = buffer.data(fph0 + 626);
    const auto *fph0_627 = buffer.data(fph0 + 627);
    const auto *fph0_629 = buffer.data(fph0 + 629);

    const auto *fpg_255 = buffer.data(fpg + 255);
    const auto *fpg_271 = buffer.data(fpg + 271);
    const auto *fpg_273 = buffer.data(fpg + 273);
    const auto *fpg_274 = buffer.data(fpg + 274);
    const auto *fpg_280 = buffer.data(fpg + 280);
    const auto *fpg_284 = buffer.data(fpg + 284);
    const auto *fpg_295 = buffer.data(fpg + 295);
    const auto *fpg_299 = buffer.data(fpg + 299);
    const auto *fpg_314 = buffer.data(fpg + 314);
    const auto *fpg_329 = buffer.data(fpg + 329);
    const auto *fpg_426 = buffer.data(fpg + 426);
    const auto *fpg_427 = buffer.data(fpg + 427);
    const auto *fpg_430 = buffer.data(fpg + 430);
    const auto *fpg_431 = buffer.data(fpg + 431);
    const auto *fpg_432 = buffer.data(fpg + 432);
    const auto *fpg_434 = buffer.data(fpg + 434);
    const auto *fpg_435 = buffer.data(fpg + 435);
    const auto *fpg_440 = buffer.data(fpg + 440);
    const auto *fpg_444 = buffer.data(fpg + 444);
    const auto *fpg_445 = buffer.data(fpg + 445);
    const auto *fpg_446 = buffer.data(fpg + 446);
    const auto *fpg_447 = buffer.data(fpg + 447);
    const auto *fpg_449 = buffer.data(fpg + 449);

    const auto *fph1_378 = buffer.data(fph1 + 378);
    const auto *fph1_379 = buffer.data(fph1 + 379);
    const auto *fph1_381 = buffer.data(fph1 + 381);
    const auto *fph1_382 = buffer.data(fph1 + 382);
    const auto *fph1_384 = buffer.data(fph1 + 384);
    const auto *fph1_385 = buffer.data(fph1 + 385);
    const auto *fph1_386 = buffer.data(fph1 + 386);
    const auto *fph1_393 = buffer.data(fph1 + 393);
    const auto *fph1_399 = buffer.data(fph1 + 399);
    const auto *fph1_400 = buffer.data(fph1 + 400);
    const auto *fph1_402 = buffer.data(fph1 + 402);
    const auto *fph1_594 = buffer.data(fph1 + 594);
    const auto *fph1_595 = buffer.data(fph1 + 595);
    const auto *fph1_603 = buffer.data(fph1 + 603);
    const auto *fph1_604 = buffer.data(fph1 + 604);
    const auto *fph1_605 = buffer.data(fph1 + 605);
    const auto *fph1_606 = buffer.data(fph1 + 606);
    const auto *fph1_608 = buffer.data(fph1 + 608);
    const auto *fph1_609 = buffer.data(fph1 + 609);
    const auto *fph1_614 = buffer.data(fph1 + 614);
    const auto *fph1_618 = buffer.data(fph1 + 618);
    const auto *fph1_624 = buffer.data(fph1 + 624);
    const auto *fph1_625 = buffer.data(fph1 + 625);
    const auto *fph1_626 = buffer.data(fph1 + 626);
    const auto *fph1_627 = buffer.data(fph1 + 627);
    const auto *fph1_629 = buffer.data(fph1 + 629);

    const auto *gsh0_194 = buffer.data(gsh0 + 194);
    const auto *gsh0_198 = buffer.data(gsh0 + 198);
    const auto *gsh0_210 = buffer.data(gsh0 + 210);
    const auto *gsh0_211 = buffer.data(gsh0 + 211);
    const auto *gsh0_213 = buffer.data(gsh0 + 213);
    const auto *gsh0_215 = buffer.data(gsh0 + 215);
    const auto *gsh0_216 = buffer.data(gsh0 + 216);
    const auto *gsh0_218 = buffer.data(gsh0 + 218);
    const auto *gsh0_219 = buffer.data(gsh0 + 219);
    const auto *gsh0_225 = buffer.data(gsh0 + 225);
    const auto *gsh0_227 = buffer.data(gsh0 + 227);
    const auto *gsh0_228 = buffer.data(gsh0 + 228);
    const auto *gsh0_230 = buffer.data(gsh0 + 230);
    const auto *gsh0_233 = buffer.data(gsh0 + 233);
    const auto *gsh0_236 = buffer.data(gsh0 + 236);
    const auto *gsh0_240 = buffer.data(gsh0 + 240);
    const auto *gsh0_248 = buffer.data(gsh0 + 248);
    const auto *gsh0_249 = buffer.data(gsh0 + 249);
    const auto *gsh0_251 = buffer.data(gsh0 + 251);

    const auto *gsg_135 = buffer.data(gsg + 135);
    const auto *gsg_140 = buffer.data(gsg + 140);
    const auto *gsg_144 = buffer.data(gsg + 144);
    const auto *gsg_149 = buffer.data(gsg + 149);
    const auto *gsg_150 = buffer.data(gsg + 150);
    const auto *gsg_151 = buffer.data(gsg + 151);
    const auto *gsg_153 = buffer.data(gsg + 153);
    const auto *gsg_155 = buffer.data(gsg + 155);
    const auto *gsg_156 = buffer.data(gsg + 156);
    const auto *gsg_158 = buffer.data(gsg + 158);
    const auto *gsg_159 = buffer.data(gsg + 159);
    const auto *gsg_160 = buffer.data(gsg + 160);
    const auto *gsg_161 = buffer.data(gsg + 161);
    const auto *gsg_162 = buffer.data(gsg + 162);
    const auto *gsg_163 = buffer.data(gsg + 163);
    const auto *gsg_164 = buffer.data(gsg + 164);
    const auto *gsg_167 = buffer.data(gsg + 167);
    const auto *gsg_170 = buffer.data(gsg + 170);
    const auto *gsg_174 = buffer.data(gsg + 174);
    const auto *gsg_175 = buffer.data(gsg + 175);
    const auto *gsg_176 = buffer.data(gsg + 176);
    const auto *gsg_177 = buffer.data(gsg + 177);
    const auto *gsg_178 = buffer.data(gsg + 178);
    const auto *gsg_179 = buffer.data(gsg + 179);

    const auto *gsh1_194 = buffer.data(gsh1 + 194);
    const auto *gsh1_198 = buffer.data(gsh1 + 198);
    const auto *gsh1_210 = buffer.data(gsh1 + 210);
    const auto *gsh1_211 = buffer.data(gsh1 + 211);
    const auto *gsh1_213 = buffer.data(gsh1 + 213);
    const auto *gsh1_215 = buffer.data(gsh1 + 215);
    const auto *gsh1_216 = buffer.data(gsh1 + 216);
    const auto *gsh1_218 = buffer.data(gsh1 + 218);
    const auto *gsh1_219 = buffer.data(gsh1 + 219);
    const auto *gsh1_225 = buffer.data(gsh1 + 225);
    const auto *gsh1_227 = buffer.data(gsh1 + 227);
    const auto *gsh1_228 = buffer.data(gsh1 + 228);
    const auto *gsh1_230 = buffer.data(gsh1 + 230);
    const auto *gsh1_233 = buffer.data(gsh1 + 233);
    const auto *gsh1_236 = buffer.data(gsh1 + 236);
    const auto *gsh1_240 = buffer.data(gsh1 + 240);
    const auto *gsh1_248 = buffer.data(gsh1 + 248);
    const auto *gsh1_249 = buffer.data(gsh1 + 249);
    const auto *gsh1_251 = buffer.data(gsh1 + 251);

    const auto *gpf0_290 = buffer.data(gpf0 + 290);
    const auto *gpf0_291 = buffer.data(gpf0 + 291);
    const auto *gpf0_292 = buffer.data(gpf0 + 292);
    const auto *gpf0_310 = buffer.data(gpf0 + 310);
    const auto *gpf0_311 = buffer.data(gpf0 + 311);
    const auto *gpf0_313 = buffer.data(gpf0 + 313);
    const auto *gpf0_315 = buffer.data(gpf0 + 315);
    const auto *gpf0_316 = buffer.data(gpf0 + 316);
    const auto *gpf0_317 = buffer.data(gpf0 + 317);
    const auto *gpf0_318 = buffer.data(gpf0 + 318);
    const auto *gpf0_319 = buffer.data(gpf0 + 319);
    const auto *gpf0_325 = buffer.data(gpf0 + 325);
    const auto *gpf0_328 = buffer.data(gpf0 + 328);
    const auto *gpf0_329 = buffer.data(gpf0 + 329);
    const auto *gpf0_342 = buffer.data(gpf0 + 342);

    const auto *gpf1_290 = buffer.data(gpf1 + 290);
    const auto *gpf1_291 = buffer.data(gpf1 + 291);
    const auto *gpf1_292 = buffer.data(gpf1 + 292);
    const auto *gpf1_310 = buffer.data(gpf1 + 310);
    const auto *gpf1_311 = buffer.data(gpf1 + 311);
    const auto *gpf1_313 = buffer.data(gpf1 + 313);
    const auto *gpf1_315 = buffer.data(gpf1 + 315);
    const auto *gpf1_316 = buffer.data(gpf1 + 316);
    const auto *gpf1_317 = buffer.data(gpf1 + 317);
    const auto *gpf1_318 = buffer.data(gpf1 + 318);
    const auto *gpf1_319 = buffer.data(gpf1 + 319);
    const auto *gpf1_325 = buffer.data(gpf1 + 325);
    const auto *gpf1_328 = buffer.data(gpf1 + 328);
    const auto *gpf1_329 = buffer.data(gpf1 + 329);
    const auto *gpf1_342 = buffer.data(gpf1 + 342);

    const auto *gpg_425 = buffer.data(gpg + 425);
    const auto *gpg_429 = buffer.data(gpg + 429);
    const auto *gpg_430 = buffer.data(gpg + 430);
    const auto *gpg_431 = buffer.data(gpg + 431);
    const auto *gpg_432 = buffer.data(gpg + 432);
    const auto *gpg_434 = buffer.data(gpg + 434);
    const auto *gpg_435 = buffer.data(gpg + 435);
    const auto *gpg_436 = buffer.data(gpg + 436);
    const auto *gpg_437 = buffer.data(gpg + 437);
    const auto *gpg_438 = buffer.data(gpg + 438);
    const auto *gpg_439 = buffer.data(gpg + 439);
    const auto *gpg_440 = buffer.data(gpg + 440);
    const auto *gpg_444 = buffer.data(gpg + 444);
    const auto *gpg_445 = buffer.data(gpg + 445);
    const auto *gpg_446 = buffer.data(gpg + 446);
    const auto *gpg_447 = buffer.data(gpg + 447);
    const auto *gpg_449 = buffer.data(gpg + 449);
    const auto *gpg_450 = buffer.data(gpg + 450);
    const auto *gpg_451 = buffer.data(gpg + 451);
    const auto *gpg_453 = buffer.data(gpg + 453);
    const auto *gpg_460 = buffer.data(gpg + 460);
    const auto *gpg_461 = buffer.data(gpg + 461);
    const auto *gpg_462 = buffer.data(gpg + 462);
    const auto *gpg_463 = buffer.data(gpg + 463);
    const auto *gpg_464 = buffer.data(gpg + 464);
    const auto *gpg_465 = buffer.data(gpg + 465);
    const auto *gpg_466 = buffer.data(gpg + 466);
    const auto *gpg_468 = buffer.data(gpg + 468);
    const auto *gpg_470 = buffer.data(gpg + 470);
    const auto *gpg_471 = buffer.data(gpg + 471);
    const auto *gpg_473 = buffer.data(gpg + 473);
    const auto *gpg_474 = buffer.data(gpg + 474);
    const auto *gpg_475 = buffer.data(gpg + 475);
    const auto *gpg_476 = buffer.data(gpg + 476);
    const auto *gpg_477 = buffer.data(gpg + 477);
    const auto *gpg_478 = buffer.data(gpg + 478);
    const auto *gpg_479 = buffer.data(gpg + 479);
    const auto *gpg_480 = buffer.data(gpg + 480);
    const auto *gpg_481 = buffer.data(gpg + 481);
    const auto *gpg_483 = buffer.data(gpg + 483);
    const auto *gpg_485 = buffer.data(gpg + 485);
    const auto *gpg_488 = buffer.data(gpg + 488);
    const auto *gpg_489 = buffer.data(gpg + 489);
    const auto *gpg_490 = buffer.data(gpg + 490);
    const auto *gpg_491 = buffer.data(gpg + 491);
    const auto *gpg_492 = buffer.data(gpg + 492);
    const auto *gpg_493 = buffer.data(gpg + 493);
    const auto *gpg_494 = buffer.data(gpg + 494);
    const auto *gpg_505 = buffer.data(gpg + 505);
    const auto *gpg_506 = buffer.data(gpg + 506);
    const auto *gpg_507 = buffer.data(gpg + 507);
    const auto *gpg_508 = buffer.data(gpg + 508);
    const auto *gpg_509 = buffer.data(gpg + 509);
    const auto *gpg_512 = buffer.data(gpg + 512);

#pragma omp simd aligned(t_593, t_594, t_595, pa_x, pb_y, pc_x, pc_y, fph0_594, fph0_595, \
                         fpg_426, fpg_427, fph1_594, fph1_595, gsh0_194, \
                         gsh1_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_593[k] = pb_y[k] * gsh0_194[k]
                   - f_9 * pc_y[k] * gsh1_194[k];

        t_594[k] = pa_x[k] * fph0_594[k]
                   + f_10 * fpg_426[k]
                   - f_9 * pc_x[k] * fph1_594[k];

        t_595[k] = pa_x[k] * fph0_595[k]
                   + f_10 * fpg_427[k]
                   - f_9 * pc_x[k] * fph1_595[k];
    }

#pragma omp simd aligned(t_596, t_597, t_598, t_599, pb_y, pc_x, pc_y, fpg_430, fpg_431, \
                         gsh0_198, gsg_140, gsh1_198, gpg_425, gpg_430, \
                         gpg_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_596[k] = f_1 * gsg_140[k]
                   + f_4 * pc_y[k] * gpg_425[k];

        t_597[k] = pb_y[k] * gsh0_198[k]
                   - f_9 * pc_y[k] * gsh1_198[k];

        t_598[k] = f_1 * fpg_430[k]
                   + f_4 * pc_x[k] * gpg_430[k];

        t_599[k] = f_1 * fpg_431[k]
                   + f_4 * pc_x[k] * gpg_431[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, t_603, pa_x, pc_x, pc_y, fph0_603, fpg_432, \
                         fpg_434, fph1_603, gsg_144, gpg_429, gpg_432, \
                         gpg_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = f_1 * fpg_432[k]
                   + f_4 * pc_x[k] * gpg_432[k];

        t_601[k] = f_1 * gsg_144[k]
                   + f_4 * pc_y[k] * gpg_429[k];

        t_602[k] = f_1 * fpg_434[k]
                   + f_4 * pc_x[k] * gpg_434[k];

        t_603[k] = pa_x[k] * fph0_603[k]
                   - f_9 * pc_x[k] * fph1_603[k];
    }

#pragma omp simd aligned(t_604, t_605, t_606, t_607, pa_x, pc_x, pc_y, fph0_604, fph0_605, \
                         fph0_606, fph1_604, fph1_605, fph1_606, gsg_149, \
                         gpg_434 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_604[k] = pa_x[k] * fph0_604[k]
                   - f_9 * pc_x[k] * fph1_604[k];

        t_605[k] = pa_x[k] * fph0_605[k]
                   - f_9 * pc_x[k] * fph1_605[k];

        t_606[k] = pa_x[k] * fph0_606[k]
                   - f_9 * pc_x[k] * fph1_606[k];

        t_607[k] = f_1 * gsg_149[k]
                   + f_4 * pc_y[k] * gpg_434[k];
    }

#pragma omp simd aligned(t_608, t_609, t_610, t_611, pa_x, pc_x, pc_y, pc_z, fph0_608, \
                         fph0_609, fpg_255, fpg_435, fph1_608, fph1_609, gsg_135, \
                         gpg_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_608[k] = pa_x[k] * fph0_608[k]
                   - f_9 * pc_x[k] * fph1_608[k];

        t_609[k] = pa_x[k] * fph0_609[k]
                   + f_12 * fpg_435[k]
                   - f_9 * pc_x[k] * fph1_609[k];

        t_610[k] = f_4 * pc_y[k] * gpg_435[k];

        t_611[k] = f_11 * fpg_255[k]
                   + f_1 * gsg_135[k]
                   + f_4 * pc_z[k] * gpg_435[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, pa_x, pc_x, pc_y, fph0_614, fpg_440, fph1_614, \
                         gpf0_290, gpf1_290, gpg_436, gpg_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = f_5 * gpf0_290[k]
                   - f_6 * gpf1_290[k]
                   + f_4 * pc_y[k] * gpg_436[k];

        t_613[k] = f_4 * pc_y[k] * gpg_437[k];

        t_614[k] = pa_x[k] * fph0_614[k]
                   + f_11 * fpg_440[k]
                   - f_9 * pc_x[k] * fph1_614[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, pc_y, gpf0_291, gpf0_292, gpf1_291, gpf1_292, \
                         gpg_438, gpg_439, gpg_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = f_7 * gpf0_291[k]
                   - f_8 * gpf1_291[k]
                   + f_4 * pc_y[k] * gpg_438[k];

        t_616[k] = f_5 * gpf0_292[k]
                   - f_6 * gpf1_292[k]
                   + f_4 * pc_y[k] * gpg_439[k];

        t_617[k] = f_4 * pc_y[k] * gpg_440[k];
    }

#pragma omp simd aligned(t_618, t_619, t_620, t_621, pa_x, pc_x, fph0_618, fpg_444, fpg_445, \
                         fpg_446, fpg_447, fph1_618, gpg_445, gpg_446, \
                         gpg_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_618[k] = pa_x[k] * fph0_618[k]
                   + f_10 * fpg_444[k]
                   - f_9 * pc_x[k] * fph1_618[k];

        t_619[k] = f_1 * fpg_445[k]
                   + f_4 * pc_x[k] * gpg_445[k];

        t_620[k] = f_1 * fpg_446[k]
                   + f_4 * pc_x[k] * gpg_446[k];

        t_621[k] = f_1 * fpg_447[k]
                   + f_4 * pc_x[k] * gpg_447[k];
    }

#pragma omp simd aligned(t_622, t_623, t_624, t_625, pa_x, pc_x, pc_y, fph0_624, fph0_625, \
                         fpg_449, fph1_624, fph1_625, gpg_444, \
                         gpg_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_622[k] = f_4 * pc_y[k] * gpg_444[k];

        t_623[k] = f_1 * fpg_449[k]
                   + f_4 * pc_x[k] * gpg_449[k];

        t_624[k] = pa_x[k] * fph0_624[k]
                   - f_9 * pc_x[k] * fph1_624[k];

        t_625[k] = pa_x[k] * fph0_625[k]
                   - f_9 * pc_x[k] * fph1_625[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, t_629, pa_x, pc_x, pc_y, fph0_626, fph0_627, \
                         fph0_629, fph1_626, fph1_627, fph1_629, \
                         gpg_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = pa_x[k] * fph0_626[k]
                   - f_9 * pc_x[k] * fph1_626[k];

        t_627[k] = pa_x[k] * fph0_627[k]
                   - f_9 * pc_x[k] * fph1_627[k];

        t_628[k] = f_4 * pc_y[k] * gpg_449[k];

        t_629[k] = pa_x[k] * fph0_629[k]
                   - f_9 * pc_x[k] * fph1_629[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, pb_x, pc_x, pc_z, gsh0_210, gsh0_211, gsg_150, \
                         gsg_151, gsh1_210, gsh1_211, gpg_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = pb_x[k] * gsh0_210[k]
                   + f_12 * gsg_150[k]
                   - f_9 * pc_x[k] * gsh1_210[k];

        t_631[k] = pb_x[k] * gsh0_211[k]
                   + f_0 * gsg_151[k]
                   - f_9 * pc_x[k] * gsh1_211[k];

        t_632[k] = f_4 * pc_z[k] * gpg_450[k];
    }

#pragma omp simd aligned(t_633, t_634, t_635, pb_x, pc_x, pc_z, gsh0_213, gsh0_215, gsg_153, \
                         gsg_155, gsh1_213, gsh1_215, gpg_451 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_633[k] = pb_x[k] * gsh0_213[k]
                   + f_11 * gsg_153[k]
                   - f_9 * pc_x[k] * gsh1_213[k];

        t_634[k] = f_4 * pc_z[k] * gpg_451[k];

        t_635[k] = pb_x[k] * gsh0_215[k]
                   + f_11 * gsg_155[k]
                   - f_9 * pc_x[k] * gsh1_215[k];
    }

#pragma omp simd aligned(t_636, t_637, t_638, pb_x, pc_x, pc_z, gsh0_216, gsh0_218, gsg_156, \
                         gsg_158, gsh1_216, gsh1_218, gpg_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_636[k] = pb_x[k] * gsh0_216[k]
                   + f_10 * gsg_156[k]
                   - f_9 * pc_x[k] * gsh1_216[k];

        t_637[k] = f_4 * pc_z[k] * gpg_453[k];

        t_638[k] = pb_x[k] * gsh0_218[k]
                   + f_10 * gsg_158[k]
                   - f_9 * pc_x[k] * gsh1_218[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, t_642, pb_x, pc_x, gsh0_219, gsg_159, gsg_160, \
                         gsg_161, gsg_162, gsh1_219, gpg_460, gpg_461, \
                         gpg_462 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = pb_x[k] * gsh0_219[k]
                   + f_10 * gsg_159[k]
                   - f_9 * pc_x[k] * gsh1_219[k];

        t_640[k] = f_1 * gsg_160[k]
                   + f_4 * pc_x[k] * gpg_460[k];

        t_641[k] = f_1 * gsg_161[k]
                   + f_4 * pc_x[k] * gpg_461[k];

        t_642[k] = f_1 * gsg_162[k]
                   + f_4 * pc_x[k] * gpg_462[k];
    }

#pragma omp simd aligned(t_643, t_644, t_645, t_646, pb_x, pc_x, pc_z, gsh0_225, gsg_163, \
                         gsg_164, gsh1_225, gpg_460, gpg_463, gpg_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_643[k] = f_1 * gsg_163[k]
                   + f_4 * pc_x[k] * gpg_463[k];

        t_644[k] = f_1 * gsg_164[k]
                   + f_4 * pc_x[k] * gpg_464[k];

        t_645[k] = pb_x[k] * gsh0_225[k]
                   - f_9 * pc_x[k] * gsh1_225[k];

        t_646[k] = f_4 * pc_z[k] * gpg_460[k];
    }

#pragma omp simd aligned(t_647, t_648, t_649, t_650, pb_x, pc_x, pc_y, fpg_284, gsh0_227, \
                         gsh0_228, gsh0_230, gsh1_227, gsh1_228, gsh1_230, \
                         gpg_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_647[k] = pb_x[k] * gsh0_227[k]
                   - f_9 * pc_x[k] * gsh1_227[k];

        t_648[k] = pb_x[k] * gsh0_228[k]
                   - f_9 * pc_x[k] * gsh1_228[k];

        t_649[k] = f_0 * fpg_284[k]
                   + f_4 * pc_y[k] * gpg_464[k];

        t_650[k] = pb_x[k] * gsh0_230[k]
                   - f_9 * pc_x[k] * gsh1_230[k];
    }

#pragma omp simd aligned(t_651, t_652, t_653, t_654, t_655, pc_x, pc_z, gpf0_310, gpf0_311, \
                         gpf0_313, gpf1_310, gpf1_311, gpf1_313, gpg_465, gpg_466, \
                         gpg_468 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_651[k] = f_2 * gpf0_310[k]
                   - f_3 * gpf1_310[k]
                   + f_4 * pc_x[k] * gpg_465[k];

        t_652[k] = f_15 * gpf0_311[k]
                   - f_16 * gpf1_311[k]
                   + f_4 * pc_x[k] * gpg_466[k];

        t_653[k] = f_4 * pc_z[k] * gpg_465[k];

        t_654[k] = f_7 * gpf0_313[k]
                   - f_8 * gpf1_313[k]
                   + f_4 * pc_x[k] * gpg_468[k];

        t_655[k] = f_4 * pc_z[k] * gpg_466[k];
    }

#pragma omp simd aligned(t_656, t_657, t_658, t_659, pc_x, pc_z, gpf0_315, gpf0_316, gpf0_318, \
                         gpf1_315, gpf1_316, gpf1_318, gpg_468, gpg_470, gpg_471, \
                         gpg_473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_656[k] = f_7 * gpf0_315[k]
                   - f_8 * gpf1_315[k]
                   + f_4 * pc_x[k] * gpg_470[k];

        t_657[k] = f_5 * gpf0_316[k]
                   - f_6 * gpf1_316[k]
                   + f_4 * pc_x[k] * gpg_471[k];

        t_658[k] = f_4 * pc_z[k] * gpg_468[k];

        t_659[k] = f_5 * gpf0_318[k]
                   - f_6 * gpf1_318[k]
                   + f_4 * pc_x[k] * gpg_473[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, t_664, t_665, pc_x, gpf0_319, gpf1_319, \
                         gpg_474, gpg_475, gpg_476, gpg_477, gpg_478, \
                         gpg_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = f_5 * gpf0_319[k]
                   - f_6 * gpf1_319[k]
                   + f_4 * pc_x[k] * gpg_474[k];

        t_661[k] = f_4 * pc_x[k] * gpg_475[k];

        t_662[k] = f_4 * pc_x[k] * gpg_476[k];

        t_663[k] = f_4 * pc_x[k] * gpg_477[k];

        t_664[k] = f_4 * pc_x[k] * gpg_478[k];

        t_665[k] = f_4 * pc_x[k] * gpg_479[k];
    }

#pragma omp simd aligned(t_666, t_667, t_668, t_669, pc_y, pc_z, fpg_295, gsg_160, gpf0_316, \
                         gpf0_317, gpf1_316, gpf1_317, gpg_475, gpg_476, \
                         gpg_477 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_666[k] = f_0 * fpg_295[k]
                   + f_1 * gsg_160[k]
                   + f_2 * gpf0_316[k]
                   - f_3 * gpf1_316[k]
                   + f_4 * pc_y[k] * gpg_475[k];

        t_667[k] = f_4 * pc_z[k] * gpg_475[k];

        t_668[k] = f_5 * gpf0_316[k]
                   - f_6 * gpf1_316[k]
                   + f_4 * pc_z[k] * gpg_476[k];

        t_669[k] = f_7 * gpf0_317[k]
                   - f_8 * gpf1_317[k]
                   + f_4 * pc_z[k] * gpg_477[k];
    }

#pragma omp simd aligned(t_670, t_671, t_672, t_673, pb_z, pc_y, pc_z, fpg_299, gsh0_210, \
                         gsh0_211, gsg_164, gsh1_210, gsh1_211, gpf0_319, gpf1_319, \
                         gpg_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_670[k] = f_0 * fpg_299[k]
                   + f_1 * gsg_164[k]
                   + f_4 * pc_y[k] * gpg_479[k];

        t_671[k] = f_2 * gpf0_319[k]
                   - f_3 * gpf1_319[k]
                   + f_4 * pc_z[k] * gpg_479[k];

        t_672[k] = pb_z[k] * gsh0_210[k]
                   - f_9 * pc_z[k] * gsh1_210[k];

        t_673[k] = pb_z[k] * gsh0_211[k]
                   - f_9 * pc_z[k] * gsh1_211[k];
    }

#pragma omp simd aligned(t_674, t_675, t_676, t_677, pb_z, pc_x, pc_z, gsh0_213, gsg_150, \
                         gsg_151, gsh1_213, gpf0_325, gpf1_325, gpg_480, gpg_481, \
                         gpg_485 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_674[k] = f_1 * gsg_150[k]
                   + f_4 * pc_z[k] * gpg_480[k];

        t_675[k] = pb_z[k] * gsh0_213[k]
                   - f_9 * pc_z[k] * gsh1_213[k];

        t_676[k] = f_1 * gsg_151[k]
                   + f_4 * pc_z[k] * gpg_481[k];

        t_677[k] = f_7 * gpf0_325[k]
                   - f_8 * gpf1_325[k]
                   + f_4 * pc_x[k] * gpg_485[k];
    }

#pragma omp simd aligned(t_678, t_679, t_680, pb_z, pc_x, pc_z, gsh0_216, gsg_153, gsh1_216, \
                         gpf0_328, gpf1_328, gpg_483, gpg_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_678[k] = pb_z[k] * gsh0_216[k]
                   - f_9 * pc_z[k] * gsh1_216[k];

        t_679[k] = f_1 * gsg_153[k]
                   + f_4 * pc_z[k] * gpg_483[k];

        t_680[k] = f_5 * gpf0_328[k]
                   - f_6 * gpf1_328[k]
                   + f_4 * pc_x[k] * gpg_488[k];
    }

#pragma omp simd aligned(t_681, t_682, t_683, t_684, t_685, t_686, pc_x, gpf0_329, gpf1_329, \
                         gpg_489, gpg_490, gpg_491, gpg_492, gpg_493, \
                         gpg_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_681[k] = f_5 * gpf0_329[k]
                   - f_6 * gpf1_329[k]
                   + f_4 * pc_x[k] * gpg_489[k];

        t_682[k] = f_4 * pc_x[k] * gpg_490[k];

        t_683[k] = f_4 * pc_x[k] * gpg_491[k];

        t_684[k] = f_4 * pc_x[k] * gpg_492[k];

        t_685[k] = f_4 * pc_x[k] * gpg_493[k];

        t_686[k] = f_4 * pc_x[k] * gpg_494[k];
    }

#pragma omp simd aligned(t_687, t_688, t_689, t_690, pb_z, pc_z, gsh0_225, gsh0_227, gsh0_228, \
                         gsg_160, gsg_161, gsg_162, gsh1_225, gsh1_227, gsh1_228, \
                         gpg_490 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_687[k] = pb_z[k] * gsh0_225[k]
                   - f_9 * pc_z[k] * gsh1_225[k];

        t_688[k] = f_1 * gsg_160[k]
                   + f_4 * pc_z[k] * gpg_490[k];

        t_689[k] = pb_z[k] * gsh0_227[k]
                   + f_10 * gsg_161[k]
                   - f_9 * pc_z[k] * gsh1_227[k];

        t_690[k] = pb_z[k] * gsh0_228[k]
                   + f_11 * gsg_162[k]
                   - f_9 * pc_z[k] * gsh1_228[k];
    }

#pragma omp simd aligned(t_691, t_692, t_693, pa_z, pb_z, pc_y, pc_z, fph0_378, fpg_314, \
                         fph1_378, gsh0_230, gsg_164, gsh1_230, \
                         gpg_494 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_691[k] = f_0 * fpg_314[k]
                   + f_4 * pc_y[k] * gpg_494[k];

        t_692[k] = pb_z[k] * gsh0_230[k]
                   + f_12 * gsg_164[k]
                   - f_9 * pc_z[k] * gsh1_230[k];

        t_693[k] = pa_z[k] * fph0_378[k]
                   - f_9 * pc_z[k] * fph1_378[k];
    }

#pragma omp simd aligned(t_694, t_695, t_696, pa_z, pb_x, pc_x, pc_z, fph0_379, fph0_381, \
                         fph1_379, fph1_381, gsh0_233, gsg_167, \
                         gsh1_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_694[k] = pa_z[k] * fph0_379[k]
                   - f_9 * pc_z[k] * fph1_379[k];

        t_695[k] = pb_x[k] * gsh0_233[k]
                   + f_0 * gsg_167[k]
                   - f_9 * pc_x[k] * gsh1_233[k];

        t_696[k] = pa_z[k] * fph0_381[k]
                   - f_9 * pc_z[k] * fph1_381[k];
    }

#pragma omp simd aligned(t_697, t_698, t_699, pa_z, pb_x, pc_x, pc_z, fph0_382, fph0_384, \
                         fpg_271, fph1_382, fph1_384, gsh0_236, gsg_170, \
                         gsh1_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_697[k] = pa_z[k] * fph0_382[k]
                   + f_1 * fpg_271[k]
                   - f_9 * pc_z[k] * fph1_382[k];

        t_698[k] = pb_x[k] * gsh0_236[k]
                   + f_11 * gsg_170[k]
                   - f_9 * pc_x[k] * gsh1_236[k];

        t_699[k] = pa_z[k] * fph0_384[k]
                   - f_9 * pc_z[k] * fph1_384[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, pa_z, pb_x, pc_x, pc_z, fph0_385, fph0_386, \
                         fpg_273, fpg_274, fph1_385, fph1_386, gsh0_240, gsg_174, \
                         gsh1_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = pa_z[k] * fph0_385[k]
                   + f_1 * fpg_273[k]
                   - f_9 * pc_z[k] * fph1_385[k];

        t_701[k] = pa_z[k] * fph0_386[k]
                   + f_10 * fpg_274[k]
                   - f_9 * pc_z[k] * fph1_386[k];

        t_702[k] = pb_x[k] * gsh0_240[k]
                   + f_10 * gsg_174[k]
                   - f_9 * pc_x[k] * gsh1_240[k];
    }

#pragma omp simd aligned(t_703, t_704, t_705, t_706, t_707, pc_x, gsg_175, gsg_176, gsg_177, \
                         gsg_178, gsg_179, gpg_505, gpg_506, gpg_507, gpg_508, \
                         gpg_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_703[k] = f_1 * gsg_175[k]
                   + f_4 * pc_x[k] * gpg_505[k];

        t_704[k] = f_1 * gsg_176[k]
                   + f_4 * pc_x[k] * gpg_506[k];

        t_705[k] = f_1 * gsg_177[k]
                   + f_4 * pc_x[k] * gpg_507[k];

        t_706[k] = f_1 * gsg_178[k]
                   + f_4 * pc_x[k] * gpg_508[k];

        t_707[k] = f_1 * gsg_179[k]
                   + f_4 * pc_x[k] * gpg_509[k];
    }

#pragma omp simd aligned(t_708, t_709, t_710, t_711, pa_z, pb_x, pc_x, pc_z, fph0_393, \
                         fpg_280, fph1_393, gsh0_248, gsh0_249, gsh1_248, gsh1_249, \
                         gpg_505 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_708[k] = pa_z[k] * fph0_393[k]
                   - f_9 * pc_z[k] * fph1_393[k];

        t_709[k] = f_1 * fpg_280[k]
                   + f_4 * pc_z[k] * gpg_505[k];

        t_710[k] = pb_x[k] * gsh0_248[k]
                   - f_9 * pc_x[k] * gsh1_248[k];

        t_711[k] = pb_x[k] * gsh0_249[k]
                   - f_9 * pc_x[k] * gsh1_249[k];
    }

#pragma omp simd aligned(t_712, t_713, t_714, pa_z, pb_x, pc_x, pc_y, pc_z, fph0_399, fpg_329, \
                         fph1_399, gsh0_251, gsh1_251, gpg_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = f_11 * fpg_329[k]
                   + f_4 * pc_y[k] * gpg_509[k];

        t_713[k] = pb_x[k] * gsh0_251[k]
                   - f_9 * pc_x[k] * gsh1_251[k];

        t_714[k] = pa_z[k] * fph0_399[k]
                   - f_9 * pc_z[k] * fph1_399[k];
    }

#pragma omp simd aligned(t_715, t_716, t_717, pa_z, pc_x, pc_z, fph0_400, fph0_402, fph1_400, \
                         fph1_402, gpf0_342, gpf1_342, gpg_512 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_715[k] = pa_z[k] * fph0_400[k]
                   - f_9 * pc_z[k] * fph1_400[k];

        t_716[k] = f_15 * gpf0_342[k]
                   - f_16 * gpf1_342[k]
                   + f_4 * pc_x[k] * gpg_512[k];

        t_717[k] = pa_z[k] * fph0_402[k]
                   - f_9 * pc_z[k] * fph1_402[k];
    }
}

static auto
compute_prim_gph_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t dph0, const size_t dph1,
                                                          const size_t fph0, const size_t fpg,
                                                          const size_t fph1, const size_t gsh0,
                                                          const size_t gsg, const size_t gsh1,
                                                          const size_t gpf0, const size_t gpf1,
                                                          const size_t gpg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 2.0 / gamma;
    const auto f_3 = 2.0 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = gamma / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 2.5 / q;
    const auto f_13 = 1.0 / p;
    const auto f_14 = gamma / (p * q);
    const auto f_15 = 1.5 / gamma;
    const auto f_16 = 1.5 * p / (gamma * q);
    const auto f_17 = 0.5 / p;
    const auto f_18 = 0.5 * gamma / (p * q);

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dph0_225 = buffer.data(dph0 + 225);
    const auto *dph0_314 = buffer.data(dph0 + 314);
    const auto *dph0_377 = buffer.data(dph0 + 377);

    const auto *dph1_225 = buffer.data(dph1 + 225);
    const auto *dph1_314 = buffer.data(dph1 + 314);
    const auto *dph1_377 = buffer.data(dph1 + 377);

    const auto *fph0_405 = buffer.data(fph0 + 405);
    const auto *fph0_414 = buffer.data(fph0 + 414);
    const auto *fph0_416 = buffer.data(fph0 + 416);
    const auto *fph0_417 = buffer.data(fph0 + 417);
    const auto *fph0_477 = buffer.data(fph0 + 477);
    const auto *fph0_503 = buffer.data(fph0 + 503);
    const auto *fph0_566 = buffer.data(fph0 + 566);
    const auto *fph0_567 = buffer.data(fph0 + 567);
    const auto *fph0_568 = buffer.data(fph0 + 568);
    const auto *fph0_569 = buffer.data(fph0 + 569);
    const auto *fph0_570 = buffer.data(fph0 + 570);
    const auto *fph0_571 = buffer.data(fph0 + 571);
    const auto *fph0_572 = buffer.data(fph0 + 572);
    const auto *fph0_573 = buffer.data(fph0 + 573);
    const auto *fph0_574 = buffer.data(fph0 + 574);
    const auto *fph0_575 = buffer.data(fph0 + 575);
    const auto *fph0_576 = buffer.data(fph0 + 576);

    const auto *fpg_295 = buffer.data(fpg + 295);
    const auto *fpg_296 = buffer.data(fpg + 296);
    const auto *fpg_297 = buffer.data(fpg + 297);
    const auto *fpg_299 = buffer.data(fpg + 299);
    const auto *fpg_310 = buffer.data(fpg + 310);
    const auto *fpg_325 = buffer.data(fpg + 325);
    const auto *fpg_340 = buffer.data(fpg + 340);
    const auto *fpg_344 = buffer.data(fpg + 344);
    const auto *fpg_355 = buffer.data(fpg + 355);
    const auto *fpg_357 = buffer.data(fpg + 357);
    const auto *fpg_358 = buffer.data(fpg + 358);
    const auto *fpg_359 = buffer.data(fpg + 359);
    const auto *fpg_374 = buffer.data(fpg + 374);
    const auto *fpg_387 = buffer.data(fpg + 387);
    const auto *fpg_388 = buffer.data(fpg + 388);
    const auto *fpg_389 = buffer.data(fpg + 389);
    const auto *fpg_400 = buffer.data(fpg + 400);
    const auto *fpg_402 = buffer.data(fpg + 402);
    const auto *fpg_403 = buffer.data(fpg + 403);
    const auto *fpg_404 = buffer.data(fpg + 404);
    const auto *fpg_405 = buffer.data(fpg + 405);
    const auto *fpg_406 = buffer.data(fpg + 406);
    const auto *fpg_407 = buffer.data(fpg + 407);
    const auto *fpg_408 = buffer.data(fpg + 408);
    const auto *fpg_409 = buffer.data(fpg + 409);
    const auto *fpg_410 = buffer.data(fpg + 410);

    const auto *fph1_405 = buffer.data(fph1 + 405);
    const auto *fph1_414 = buffer.data(fph1 + 414);
    const auto *fph1_416 = buffer.data(fph1 + 416);
    const auto *fph1_417 = buffer.data(fph1 + 417);
    const auto *fph1_477 = buffer.data(fph1 + 477);
    const auto *fph1_503 = buffer.data(fph1 + 503);
    const auto *fph1_566 = buffer.data(fph1 + 566);
    const auto *fph1_567 = buffer.data(fph1 + 567);
    const auto *fph1_568 = buffer.data(fph1 + 568);
    const auto *fph1_569 = buffer.data(fph1 + 569);
    const auto *fph1_570 = buffer.data(fph1 + 570);
    const auto *fph1_571 = buffer.data(fph1 + 571);
    const auto *fph1_572 = buffer.data(fph1 + 572);
    const auto *fph1_573 = buffer.data(fph1 + 573);
    const auto *fph1_574 = buffer.data(fph1 + 574);
    const auto *fph1_575 = buffer.data(fph1 + 575);
    const auto *fph1_576 = buffer.data(fph1 + 576);

    const auto *gsh0_252 = buffer.data(gsh0 + 252);
    const auto *gsh0_253 = buffer.data(gsh0 + 253);
    const auto *gsh0_254 = buffer.data(gsh0 + 254);
    const auto *gsh0_255 = buffer.data(gsh0 + 255);
    const auto *gsh0_256 = buffer.data(gsh0 + 256);
    const auto *gsh0_257 = buffer.data(gsh0 + 257);
    const auto *gsh0_258 = buffer.data(gsh0 + 258);
    const auto *gsh0_259 = buffer.data(gsh0 + 259);
    const auto *gsh0_260 = buffer.data(gsh0 + 260);
    const auto *gsh0_261 = buffer.data(gsh0 + 261);
    const auto *gsh0_267 = buffer.data(gsh0 + 267);
    const auto *gsh0_269 = buffer.data(gsh0 + 269);
    const auto *gsh0_270 = buffer.data(gsh0 + 270);
    const auto *gsh0_272 = buffer.data(gsh0 + 272);

    const auto *gsg_175 = buffer.data(gsg + 175);
    const auto *gsg_179 = buffer.data(gsg + 179);
    const auto *gsg_180 = buffer.data(gsg + 180);
    const auto *gsg_181 = buffer.data(gsg + 181);
    const auto *gsg_182 = buffer.data(gsg + 182);
    const auto *gsg_183 = buffer.data(gsg + 183);
    const auto *gsg_184 = buffer.data(gsg + 184);
    const auto *gsg_185 = buffer.data(gsg + 185);
    const auto *gsg_186 = buffer.data(gsg + 186);
    const auto *gsg_187 = buffer.data(gsg + 187);
    const auto *gsg_188 = buffer.data(gsg + 188);
    const auto *gsg_189 = buffer.data(gsg + 189);
    const auto *gsg_190 = buffer.data(gsg + 190);
    const auto *gsg_191 = buffer.data(gsg + 191);
    const auto *gsg_192 = buffer.data(gsg + 192);
    const auto *gsg_193 = buffer.data(gsg + 193);
    const auto *gsg_194 = buffer.data(gsg + 194);
    const auto *gsg_205 = buffer.data(gsg + 205);
    const auto *gsg_206 = buffer.data(gsg + 206);
    const auto *gsg_207 = buffer.data(gsg + 207);
    const auto *gsg_208 = buffer.data(gsg + 208);
    const auto *gsg_209 = buffer.data(gsg + 209);

    const auto *gsh1_252 = buffer.data(gsh1 + 252);
    const auto *gsh1_253 = buffer.data(gsh1 + 253);
    const auto *gsh1_254 = buffer.data(gsh1 + 254);
    const auto *gsh1_255 = buffer.data(gsh1 + 255);
    const auto *gsh1_256 = buffer.data(gsh1 + 256);
    const auto *gsh1_257 = buffer.data(gsh1 + 257);
    const auto *gsh1_258 = buffer.data(gsh1 + 258);
    const auto *gsh1_259 = buffer.data(gsh1 + 259);
    const auto *gsh1_260 = buffer.data(gsh1 + 260);
    const auto *gsh1_261 = buffer.data(gsh1 + 261);
    const auto *gsh1_267 = buffer.data(gsh1 + 267);
    const auto *gsh1_269 = buffer.data(gsh1 + 269);
    const auto *gsh1_270 = buffer.data(gsh1 + 270);
    const auto *gsh1_272 = buffer.data(gsh1 + 272);

    const auto *gpf0_344 = buffer.data(gpf0 + 344);
    const auto *gpf0_345 = buffer.data(gpf0 + 345);
    const auto *gpf0_347 = buffer.data(gpf0 + 347);
    const auto *gpf0_348 = buffer.data(gpf0 + 348);
    const auto *gpf0_349 = buffer.data(gpf0 + 349);
    const auto *gpf0_350 = buffer.data(gpf0 + 350);
    const auto *gpf0_351 = buffer.data(gpf0 + 351);
    const auto *gpf0_352 = buffer.data(gpf0 + 352);
    const auto *gpf0_353 = buffer.data(gpf0 + 353);
    const auto *gpf0_354 = buffer.data(gpf0 + 354);
    const auto *gpf0_355 = buffer.data(gpf0 + 355);
    const auto *gpf0_356 = buffer.data(gpf0 + 356);
    const auto *gpf0_357 = buffer.data(gpf0 + 357);
    const auto *gpf0_358 = buffer.data(gpf0 + 358);
    const auto *gpf0_359 = buffer.data(gpf0 + 359);
    const auto *gpf0_370 = buffer.data(gpf0 + 370);
    const auto *gpf0_371 = buffer.data(gpf0 + 371);
    const auto *gpf0_372 = buffer.data(gpf0 + 372);
    const auto *gpf0_373 = buffer.data(gpf0 + 373);
    const auto *gpf0_374 = buffer.data(gpf0 + 374);
    const auto *gpf0_375 = buffer.data(gpf0 + 375);
    const auto *gpf0_376 = buffer.data(gpf0 + 376);
    const auto *gpf0_377 = buffer.data(gpf0 + 377);
    const auto *gpf0_378 = buffer.data(gpf0 + 378);
    const auto *gpf0_379 = buffer.data(gpf0 + 379);
    const auto *gpf0_380 = buffer.data(gpf0 + 380);
    const auto *gpf0_381 = buffer.data(gpf0 + 381);
    const auto *gpf0_382 = buffer.data(gpf0 + 382);
    const auto *gpf0_383 = buffer.data(gpf0 + 383);
    const auto *gpf0_384 = buffer.data(gpf0 + 384);
    const auto *gpf0_385 = buffer.data(gpf0 + 385);
    const auto *gpf0_386 = buffer.data(gpf0 + 386);
    const auto *gpf0_387 = buffer.data(gpf0 + 387);
    const auto *gpf0_388 = buffer.data(gpf0 + 388);
    const auto *gpf0_389 = buffer.data(gpf0 + 389);

    const auto *gpf1_344 = buffer.data(gpf1 + 344);
    const auto *gpf1_345 = buffer.data(gpf1 + 345);
    const auto *gpf1_347 = buffer.data(gpf1 + 347);
    const auto *gpf1_348 = buffer.data(gpf1 + 348);
    const auto *gpf1_349 = buffer.data(gpf1 + 349);
    const auto *gpf1_350 = buffer.data(gpf1 + 350);
    const auto *gpf1_351 = buffer.data(gpf1 + 351);
    const auto *gpf1_352 = buffer.data(gpf1 + 352);
    const auto *gpf1_353 = buffer.data(gpf1 + 353);
    const auto *gpf1_354 = buffer.data(gpf1 + 354);
    const auto *gpf1_355 = buffer.data(gpf1 + 355);
    const auto *gpf1_356 = buffer.data(gpf1 + 356);
    const auto *gpf1_357 = buffer.data(gpf1 + 357);
    const auto *gpf1_358 = buffer.data(gpf1 + 358);
    const auto *gpf1_359 = buffer.data(gpf1 + 359);
    const auto *gpf1_370 = buffer.data(gpf1 + 370);
    const auto *gpf1_371 = buffer.data(gpf1 + 371);
    const auto *gpf1_372 = buffer.data(gpf1 + 372);
    const auto *gpf1_373 = buffer.data(gpf1 + 373);
    const auto *gpf1_374 = buffer.data(gpf1 + 374);
    const auto *gpf1_375 = buffer.data(gpf1 + 375);
    const auto *gpf1_376 = buffer.data(gpf1 + 376);
    const auto *gpf1_377 = buffer.data(gpf1 + 377);
    const auto *gpf1_378 = buffer.data(gpf1 + 378);
    const auto *gpf1_379 = buffer.data(gpf1 + 379);
    const auto *gpf1_380 = buffer.data(gpf1 + 380);
    const auto *gpf1_381 = buffer.data(gpf1 + 381);
    const auto *gpf1_382 = buffer.data(gpf1 + 382);
    const auto *gpf1_383 = buffer.data(gpf1 + 383);
    const auto *gpf1_384 = buffer.data(gpf1 + 384);
    const auto *gpf1_385 = buffer.data(gpf1 + 385);
    const auto *gpf1_386 = buffer.data(gpf1 + 386);
    const auto *gpf1_387 = buffer.data(gpf1 + 387);
    const auto *gpf1_388 = buffer.data(gpf1 + 388);
    const auto *gpf1_389 = buffer.data(gpf1 + 389);

    const auto *gpg_514 = buffer.data(gpg + 514);
    const auto *gpg_515 = buffer.data(gpg + 515);
    const auto *gpg_517 = buffer.data(gpg + 517);
    const auto *gpg_518 = buffer.data(gpg + 518);
    const auto *gpg_519 = buffer.data(gpg + 519);
    const auto *gpg_520 = buffer.data(gpg + 520);
    const auto *gpg_521 = buffer.data(gpg + 521);
    const auto *gpg_522 = buffer.data(gpg + 522);
    const auto *gpg_523 = buffer.data(gpg + 523);
    const auto *gpg_524 = buffer.data(gpg + 524);
    const auto *gpg_525 = buffer.data(gpg + 525);
    const auto *gpg_526 = buffer.data(gpg + 526);
    const auto *gpg_527 = buffer.data(gpg + 527);
    const auto *gpg_528 = buffer.data(gpg + 528);
    const auto *gpg_529 = buffer.data(gpg + 529);
    const auto *gpg_530 = buffer.data(gpg + 530);
    const auto *gpg_531 = buffer.data(gpg + 531);
    const auto *gpg_532 = buffer.data(gpg + 532);
    const auto *gpg_533 = buffer.data(gpg + 533);
    const auto *gpg_534 = buffer.data(gpg + 534);
    const auto *gpg_535 = buffer.data(gpg + 535);
    const auto *gpg_536 = buffer.data(gpg + 536);
    const auto *gpg_537 = buffer.data(gpg + 537);
    const auto *gpg_538 = buffer.data(gpg + 538);
    const auto *gpg_539 = buffer.data(gpg + 539);
    const auto *gpg_550 = buffer.data(gpg + 550);
    const auto *gpg_551 = buffer.data(gpg + 551);
    const auto *gpg_552 = buffer.data(gpg + 552);
    const auto *gpg_553 = buffer.data(gpg + 553);
    const auto *gpg_554 = buffer.data(gpg + 554);
    const auto *gpg_555 = buffer.data(gpg + 555);
    const auto *gpg_556 = buffer.data(gpg + 556);
    const auto *gpg_557 = buffer.data(gpg + 557);
    const auto *gpg_558 = buffer.data(gpg + 558);
    const auto *gpg_559 = buffer.data(gpg + 559);
    const auto *gpg_560 = buffer.data(gpg + 560);
    const auto *gpg_561 = buffer.data(gpg + 561);
    const auto *gpg_562 = buffer.data(gpg + 562);
    const auto *gpg_563 = buffer.data(gpg + 563);
    const auto *gpg_564 = buffer.data(gpg + 564);
    const auto *gpg_565 = buffer.data(gpg + 565);
    const auto *gpg_566 = buffer.data(gpg + 566);
    const auto *gpg_567 = buffer.data(gpg + 567);
    const auto *gpg_568 = buffer.data(gpg + 568);
    const auto *gpg_569 = buffer.data(gpg + 569);
    const auto *gpg_570 = buffer.data(gpg + 570);
    const auto *gpg_571 = buffer.data(gpg + 571);
    const auto *gpg_572 = buffer.data(gpg + 572);
    const auto *gpg_573 = buffer.data(gpg + 573);
    const auto *gpg_574 = buffer.data(gpg + 574);
    const auto *gpg_575 = buffer.data(gpg + 575);
    const auto *gpg_576 = buffer.data(gpg + 576);
    const auto *gpg_577 = buffer.data(gpg + 577);
    const auto *gpg_578 = buffer.data(gpg + 578);
    const auto *gpg_579 = buffer.data(gpg + 579);
    const auto *gpg_580 = buffer.data(gpg + 580);
    const auto *gpg_581 = buffer.data(gpg + 581);
    const auto *gpg_582 = buffer.data(gpg + 582);
    const auto *gpg_583 = buffer.data(gpg + 583);
    const auto *gpg_584 = buffer.data(gpg + 584);
    const auto *gpg_595 = buffer.data(gpg + 595);
    const auto *gpg_596 = buffer.data(gpg + 596);
    const auto *gpg_597 = buffer.data(gpg + 597);
    const auto *gpg_598 = buffer.data(gpg + 598);
    const auto *gpg_599 = buffer.data(gpg + 599);

#pragma omp simd aligned(t_718, t_719, t_720, pa_z, pc_x, pc_z, fph0_405, fph1_405, gpf0_344, \
                         gpf0_345, gpf1_344, gpf1_345, gpg_514, \
                         gpg_515 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_718[k] = f_7 * gpf0_344[k]
                   - f_8 * gpf1_344[k]
                   + f_4 * pc_x[k] * gpg_514[k];

        t_719[k] = f_7 * gpf0_345[k]
                   - f_8 * gpf1_345[k]
                   + f_4 * pc_x[k] * gpg_515[k];

        t_720[k] = pa_z[k] * fph0_405[k]
                   - f_9 * pc_z[k] * fph1_405[k];
    }

#pragma omp simd aligned(t_721, t_722, t_723, t_724, pc_x, gpf0_347, gpf0_348, gpf0_349, \
                         gpf1_347, gpf1_348, gpf1_349, gpg_517, gpg_518, gpg_519, \
                         gpg_520 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_721[k] = f_5 * gpf0_347[k]
                   - f_6 * gpf1_347[k]
                   + f_4 * pc_x[k] * gpg_517[k];

        t_722[k] = f_5 * gpf0_348[k]
                   - f_6 * gpf1_348[k]
                   + f_4 * pc_x[k] * gpg_518[k];

        t_723[k] = f_5 * gpf0_349[k]
                   - f_6 * gpf1_349[k]
                   + f_4 * pc_x[k] * gpg_519[k];

        t_724[k] = f_4 * pc_x[k] * gpg_520[k];
    }

#pragma omp simd aligned(t_725, t_726, t_727, t_728, t_729, pa_z, pc_x, pc_z, fph0_414, \
                         fph1_414, gpg_521, gpg_522, gpg_523, gpg_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_725[k] = f_4 * pc_x[k] * gpg_521[k];

        t_726[k] = f_4 * pc_x[k] * gpg_522[k];

        t_727[k] = f_4 * pc_x[k] * gpg_523[k];

        t_728[k] = f_4 * pc_x[k] * gpg_524[k];

        t_729[k] = pa_z[k] * fph0_414[k]
                   - f_9 * pc_z[k] * fph1_414[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, pa_z, pc_z, fph0_416, fph0_417, fpg_295, \
                         fpg_296, fpg_297, fph1_416, fph1_417, \
                         gpg_520 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = f_1 * fpg_295[k]
                   + f_4 * pc_z[k] * gpg_520[k];

        t_731[k] = pa_z[k] * fph0_416[k]
                   + f_10 * fpg_296[k]
                   - f_9 * pc_z[k] * fph1_416[k];

        t_732[k] = pa_z[k] * fph0_417[k]
                   + f_11 * fpg_297[k]
                   - f_9 * pc_z[k] * fph1_417[k];
    }

#pragma omp simd aligned(t_733, t_734, t_735, pc_x, pc_y, pc_z, fpg_299, fpg_344, gsg_179, \
                         gpf0_349, gpf0_350, gpf1_349, gpf1_350, gpg_524, \
                         gpg_525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_733[k] = f_11 * fpg_344[k]
                   + f_1 * gsg_179[k]
                   + f_4 * pc_y[k] * gpg_524[k];

        t_734[k] = f_1 * fpg_299[k]
                   + f_2 * gpf0_349[k]
                   - f_3 * gpf1_349[k]
                   + f_4 * pc_z[k] * gpg_524[k];

        t_735[k] = f_2 * gpf0_350[k]
                   - f_3 * gpf1_350[k]
                   + f_4 * pc_x[k] * gpg_525[k];
    }

#pragma omp simd aligned(t_736, t_737, t_738, pc_x, gpf0_351, gpf0_352, gpf0_353, gpf1_351, \
                         gpf1_352, gpf1_353, gpg_526, gpg_527, \
                         gpg_528 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_736[k] = f_15 * gpf0_351[k]
                   - f_16 * gpf1_351[k]
                   + f_4 * pc_x[k] * gpg_526[k];

        t_737[k] = f_15 * gpf0_352[k]
                   - f_16 * gpf1_352[k]
                   + f_4 * pc_x[k] * gpg_527[k];

        t_738[k] = f_7 * gpf0_353[k]
                   - f_8 * gpf1_353[k]
                   + f_4 * pc_x[k] * gpg_528[k];
    }

#pragma omp simd aligned(t_739, t_740, t_741, pc_x, gpf0_354, gpf0_355, gpf0_356, gpf1_354, \
                         gpf1_355, gpf1_356, gpg_529, gpg_530, \
                         gpg_531 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_739[k] = f_7 * gpf0_354[k]
                   - f_8 * gpf1_354[k]
                   + f_4 * pc_x[k] * gpg_529[k];

        t_740[k] = f_7 * gpf0_355[k]
                   - f_8 * gpf1_355[k]
                   + f_4 * pc_x[k] * gpg_530[k];

        t_741[k] = f_5 * gpf0_356[k]
                   - f_6 * gpf1_356[k]
                   + f_4 * pc_x[k] * gpg_531[k];
    }

#pragma omp simd aligned(t_742, t_743, t_744, t_745, pc_x, gpf0_357, gpf0_358, gpf0_359, \
                         gpf1_357, gpf1_358, gpf1_359, gpg_532, gpg_533, gpg_534, \
                         gpg_535 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_742[k] = f_5 * gpf0_357[k]
                   - f_6 * gpf1_357[k]
                   + f_4 * pc_x[k] * gpg_532[k];

        t_743[k] = f_5 * gpf0_358[k]
                   - f_6 * gpf1_358[k]
                   + f_4 * pc_x[k] * gpg_533[k];

        t_744[k] = f_5 * gpf0_359[k]
                   - f_6 * gpf1_359[k]
                   + f_4 * pc_x[k] * gpg_534[k];

        t_745[k] = f_4 * pc_x[k] * gpg_535[k];
    }

#pragma omp simd aligned(t_746, t_747, t_748, t_749, t_750, pc_x, pc_y, fpg_355, gpf0_356, \
                         gpf1_356, gpg_535, gpg_536, gpg_537, gpg_538, \
                         gpg_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_746[k] = f_4 * pc_x[k] * gpg_536[k];

        t_747[k] = f_4 * pc_x[k] * gpg_537[k];

        t_748[k] = f_4 * pc_x[k] * gpg_538[k];

        t_749[k] = f_4 * pc_x[k] * gpg_539[k];

        t_750[k] = f_11 * fpg_355[k]
                   + f_2 * gpf0_356[k]
                   - f_3 * gpf1_356[k]
                   + f_4 * pc_y[k] * gpg_535[k];
    }

#pragma omp simd aligned(t_751, t_752, t_753, pc_y, pc_z, fpg_310, fpg_357, fpg_358, gsg_175, \
                         gpf0_358, gpf0_359, gpf1_358, gpf1_359, gpg_535, gpg_537, \
                         gpg_538 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_751[k] = f_1 * fpg_310[k]
                   + f_1 * gsg_175[k]
                   + f_4 * pc_z[k] * gpg_535[k];

        t_752[k] = f_11 * fpg_357[k]
                   + f_7 * gpf0_358[k]
                   - f_8 * gpf1_358[k]
                   + f_4 * pc_y[k] * gpg_537[k];

        t_753[k] = f_11 * fpg_358[k]
                   + f_5 * gpf0_359[k]
                   - f_6 * gpf1_359[k]
                   + f_4 * pc_y[k] * gpg_538[k];
    }

#pragma omp simd aligned(t_754, t_755, t_756, pa_y, pb_x, pc_x, pc_y, dph0_314, dph1_314, \
                         fph0_503, fpg_359, fph1_503, gsh0_252, gsg_180, gsh1_252, \
                         gpg_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_754[k] = f_11 * fpg_359[k]
                   + f_4 * pc_y[k] * gpg_539[k];

        t_755[k] = f_13 * dph0_314[k]
                   - f_14 * dph1_314[k]
                   + pa_y[k] * fph0_503[k]
                   - f_9 * pc_y[k] * fph1_503[k];

        t_756[k] = pb_x[k] * gsh0_252[k]
                   + f_12 * gsg_180[k]
                   - f_9 * pc_x[k] * gsh1_252[k];
    }

#pragma omp simd aligned(t_757, t_758, t_759, pb_x, pc_x, gsh0_253, gsh0_254, gsh0_255, \
                         gsg_181, gsg_182, gsg_183, gsh1_253, gsh1_254, \
                         gsh1_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_757[k] = pb_x[k] * gsh0_253[k]
                   + f_0 * gsg_181[k]
                   - f_9 * pc_x[k] * gsh1_253[k];

        t_758[k] = pb_x[k] * gsh0_254[k]
                   + f_0 * gsg_182[k]
                   - f_9 * pc_x[k] * gsh1_254[k];

        t_759[k] = pb_x[k] * gsh0_255[k]
                   + f_11 * gsg_183[k]
                   - f_9 * pc_x[k] * gsh1_255[k];
    }

#pragma omp simd aligned(t_760, t_761, t_762, pb_x, pc_x, gsh0_256, gsh0_257, gsh0_258, \
                         gsg_184, gsg_185, gsg_186, gsh1_256, gsh1_257, \
                         gsh1_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_760[k] = pb_x[k] * gsh0_256[k]
                   + f_11 * gsg_184[k]
                   - f_9 * pc_x[k] * gsh1_256[k];

        t_761[k] = pb_x[k] * gsh0_257[k]
                   + f_11 * gsg_185[k]
                   - f_9 * pc_x[k] * gsh1_257[k];

        t_762[k] = pb_x[k] * gsh0_258[k]
                   + f_10 * gsg_186[k]
                   - f_9 * pc_x[k] * gsh1_258[k];
    }

#pragma omp simd aligned(t_763, t_764, t_765, pb_x, pc_x, gsh0_259, gsh0_260, gsh0_261, \
                         gsg_187, gsg_188, gsg_189, gsh1_259, gsh1_260, \
                         gsh1_261 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_763[k] = pb_x[k] * gsh0_259[k]
                   + f_10 * gsg_187[k]
                   - f_9 * pc_x[k] * gsh1_259[k];

        t_764[k] = pb_x[k] * gsh0_260[k]
                   + f_10 * gsg_188[k]
                   - f_9 * pc_x[k] * gsh1_260[k];

        t_765[k] = pb_x[k] * gsh0_261[k]
                   + f_10 * gsg_189[k]
                   - f_9 * pc_x[k] * gsh1_261[k];
    }

#pragma omp simd aligned(t_766, t_767, t_768, t_769, t_770, pc_x, gsg_190, gsg_191, gsg_192, \
                         gsg_193, gsg_194, gpg_550, gpg_551, gpg_552, gpg_553, \
                         gpg_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_766[k] = f_1 * gsg_190[k]
                   + f_4 * pc_x[k] * gpg_550[k];

        t_767[k] = f_1 * gsg_191[k]
                   + f_4 * pc_x[k] * gpg_551[k];

        t_768[k] = f_1 * gsg_192[k]
                   + f_4 * pc_x[k] * gpg_552[k];

        t_769[k] = f_1 * gsg_193[k]
                   + f_4 * pc_x[k] * gpg_553[k];

        t_770[k] = f_1 * gsg_194[k]
                   + f_4 * pc_x[k] * gpg_554[k];
    }

#pragma omp simd aligned(t_771, t_772, t_773, t_774, pb_x, pc_x, pc_z, fpg_325, gsh0_267, \
                         gsh0_269, gsh0_270, gsh1_267, gsh1_269, gsh1_270, \
                         gpg_550 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_771[k] = pb_x[k] * gsh0_267[k]
                   - f_9 * pc_x[k] * gsh1_267[k];

        t_772[k] = f_10 * fpg_325[k]
                   + f_4 * pc_z[k] * gpg_550[k];

        t_773[k] = pb_x[k] * gsh0_269[k]
                   - f_9 * pc_x[k] * gsh1_269[k];

        t_774[k] = pb_x[k] * gsh0_270[k]
                   - f_9 * pc_x[k] * gsh1_270[k];
    }

#pragma omp simd aligned(t_775, t_776, t_777, pb_x, pc_x, pc_y, fpg_374, gsh0_272, gsh1_272, \
                         gpf0_370, gpf1_370, gpg_554, gpg_555 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_775[k] = f_10 * fpg_374[k]
                   + f_4 * pc_y[k] * gpg_554[k];

        t_776[k] = pb_x[k] * gsh0_272[k]
                   - f_9 * pc_x[k] * gsh1_272[k];

        t_777[k] = f_2 * gpf0_370[k]
                   - f_3 * gpf1_370[k]
                   + f_4 * pc_x[k] * gpg_555[k];
    }

#pragma omp simd aligned(t_778, t_779, t_780, pc_x, gpf0_371, gpf0_372, gpf0_373, gpf1_371, \
                         gpf1_372, gpf1_373, gpg_556, gpg_557, \
                         gpg_558 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_778[k] = f_15 * gpf0_371[k]
                   - f_16 * gpf1_371[k]
                   + f_4 * pc_x[k] * gpg_556[k];

        t_779[k] = f_15 * gpf0_372[k]
                   - f_16 * gpf1_372[k]
                   + f_4 * pc_x[k] * gpg_557[k];

        t_780[k] = f_7 * gpf0_373[k]
                   - f_8 * gpf1_373[k]
                   + f_4 * pc_x[k] * gpg_558[k];
    }

#pragma omp simd aligned(t_781, t_782, t_783, pc_x, gpf0_374, gpf0_375, gpf0_376, gpf1_374, \
                         gpf1_375, gpf1_376, gpg_559, gpg_560, \
                         gpg_561 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_781[k] = f_7 * gpf0_374[k]
                   - f_8 * gpf1_374[k]
                   + f_4 * pc_x[k] * gpg_559[k];

        t_782[k] = f_7 * gpf0_375[k]
                   - f_8 * gpf1_375[k]
                   + f_4 * pc_x[k] * gpg_560[k];

        t_783[k] = f_5 * gpf0_376[k]
                   - f_6 * gpf1_376[k]
                   + f_4 * pc_x[k] * gpg_561[k];
    }

#pragma omp simd aligned(t_784, t_785, t_786, t_787, pc_x, gpf0_377, gpf0_378, gpf0_379, \
                         gpf1_377, gpf1_378, gpf1_379, gpg_562, gpg_563, gpg_564, \
                         gpg_565 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_784[k] = f_5 * gpf0_377[k]
                   - f_6 * gpf1_377[k]
                   + f_4 * pc_x[k] * gpg_562[k];

        t_785[k] = f_5 * gpf0_378[k]
                   - f_6 * gpf1_378[k]
                   + f_4 * pc_x[k] * gpg_563[k];

        t_786[k] = f_5 * gpf0_379[k]
                   - f_6 * gpf1_379[k]
                   + f_4 * pc_x[k] * gpg_564[k];

        t_787[k] = f_4 * pc_x[k] * gpg_565[k];
    }

#pragma omp simd aligned(t_788, t_789, t_790, t_791, t_792, pa_z, pc_x, pc_z, dph0_225, \
                         dph1_225, fph0_477, fph1_477, gpg_566, gpg_567, gpg_568, \
                         gpg_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_788[k] = f_4 * pc_x[k] * gpg_566[k];

        t_789[k] = f_4 * pc_x[k] * gpg_567[k];

        t_790[k] = f_4 * pc_x[k] * gpg_568[k];

        t_791[k] = f_4 * pc_x[k] * gpg_569[k];

        t_792[k] = f_17 * dph0_225[k]
                   - f_18 * dph1_225[k]
                   + pa_z[k] * fph0_477[k]
                   - f_9 * pc_z[k] * fph1_477[k];
    }

#pragma omp simd aligned(t_793, t_794, pc_y, pc_z, fpg_340, fpg_387, gsg_192, gpf0_378, \
                         gpf1_378, gpg_565, gpg_567 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_793[k] = f_10 * fpg_340[k]
                   + f_4 * pc_z[k] * gpg_565[k];

        t_794[k] = f_10 * fpg_387[k]
                   + f_1 * gsg_192[k]
                   + f_7 * gpf0_378[k]
                   - f_8 * gpf1_378[k]
                   + f_4 * pc_y[k] * gpg_567[k];
    }

#pragma omp simd aligned(t_795, t_796, t_797, pc_y, pc_z, fpg_344, fpg_388, fpg_389, gsg_193, \
                         gsg_194, gpf0_379, gpf1_379, gpg_568, \
                         gpg_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_795[k] = f_10 * fpg_388[k]
                   + f_1 * gsg_193[k]
                   + f_5 * gpf0_379[k]
                   - f_6 * gpf1_379[k]
                   + f_4 * pc_y[k] * gpg_568[k];

        t_796[k] = f_10 * fpg_389[k]
                   + f_1 * gsg_194[k]
                   + f_4 * pc_y[k] * gpg_569[k];

        t_797[k] = f_10 * fpg_344[k]
                   + f_2 * gpf0_379[k]
                   - f_3 * gpf1_379[k]
                   + f_4 * pc_z[k] * gpg_569[k];
    }

#pragma omp simd aligned(t_798, t_799, t_800, pc_x, gpf0_380, gpf0_381, gpf0_382, gpf1_380, \
                         gpf1_381, gpf1_382, gpg_570, gpg_571, \
                         gpg_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_798[k] = f_2 * gpf0_380[k]
                   - f_3 * gpf1_380[k]
                   + f_4 * pc_x[k] * gpg_570[k];

        t_799[k] = f_15 * gpf0_381[k]
                   - f_16 * gpf1_381[k]
                   + f_4 * pc_x[k] * gpg_571[k];

        t_800[k] = f_15 * gpf0_382[k]
                   - f_16 * gpf1_382[k]
                   + f_4 * pc_x[k] * gpg_572[k];
    }

#pragma omp simd aligned(t_801, t_802, t_803, pc_x, gpf0_383, gpf0_384, gpf0_385, gpf1_383, \
                         gpf1_384, gpf1_385, gpg_573, gpg_574, \
                         gpg_575 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_801[k] = f_7 * gpf0_383[k]
                   - f_8 * gpf1_383[k]
                   + f_4 * pc_x[k] * gpg_573[k];

        t_802[k] = f_7 * gpf0_384[k]
                   - f_8 * gpf1_384[k]
                   + f_4 * pc_x[k] * gpg_574[k];

        t_803[k] = f_7 * gpf0_385[k]
                   - f_8 * gpf1_385[k]
                   + f_4 * pc_x[k] * gpg_575[k];
    }

#pragma omp simd aligned(t_804, t_805, t_806, pc_x, gpf0_386, gpf0_387, gpf0_388, gpf1_386, \
                         gpf1_387, gpf1_388, gpg_576, gpg_577, \
                         gpg_578 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_804[k] = f_5 * gpf0_386[k]
                   - f_6 * gpf1_386[k]
                   + f_4 * pc_x[k] * gpg_576[k];

        t_805[k] = f_5 * gpf0_387[k]
                   - f_6 * gpf1_387[k]
                   + f_4 * pc_x[k] * gpg_577[k];

        t_806[k] = f_5 * gpf0_388[k]
                   - f_6 * gpf1_388[k]
                   + f_4 * pc_x[k] * gpg_578[k];
    }

#pragma omp simd aligned(t_807, t_808, t_809, t_810, t_811, t_812, pc_x, gpf0_389, gpf1_389, \
                         gpg_579, gpg_580, gpg_581, gpg_582, gpg_583, \
                         gpg_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_807[k] = f_5 * gpf0_389[k]
                   - f_6 * gpf1_389[k]
                   + f_4 * pc_x[k] * gpg_579[k];

        t_808[k] = f_4 * pc_x[k] * gpg_580[k];

        t_809[k] = f_4 * pc_x[k] * gpg_581[k];

        t_810[k] = f_4 * pc_x[k] * gpg_582[k];

        t_811[k] = f_4 * pc_x[k] * gpg_583[k];

        t_812[k] = f_4 * pc_x[k] * gpg_584[k];
    }

#pragma omp simd aligned(t_813, t_814, t_815, pc_y, pc_z, fpg_355, fpg_400, fpg_402, gsg_190, \
                         gpf0_386, gpf0_388, gpf1_386, gpf1_388, gpg_580, \
                         gpg_582 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_813[k] = f_10 * fpg_400[k]
                   + f_2 * gpf0_386[k]
                   - f_3 * gpf1_386[k]
                   + f_4 * pc_y[k] * gpg_580[k];

        t_814[k] = f_10 * fpg_355[k]
                   + f_1 * gsg_190[k]
                   + f_4 * pc_z[k] * gpg_580[k];

        t_815[k] = f_10 * fpg_402[k]
                   + f_7 * gpf0_388[k]
                   - f_8 * gpf1_388[k]
                   + f_4 * pc_y[k] * gpg_582[k];
    }

#pragma omp simd aligned(t_816, t_817, t_818, pa_y, pc_y, dph0_377, dph1_377, fph0_566, \
                         fpg_403, fpg_404, fph1_566, gpf0_389, gpf1_389, gpg_583, \
                         gpg_584 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_816[k] = f_10 * fpg_403[k]
                   + f_5 * gpf0_389[k]
                   - f_6 * gpf1_389[k]
                   + f_4 * pc_y[k] * gpg_583[k];

        t_817[k] = f_10 * fpg_404[k]
                   + f_4 * pc_y[k] * gpg_584[k];

        t_818[k] = f_17 * dph0_377[k]
                   - f_18 * dph1_377[k]
                   + pa_y[k] * fph0_566[k]
                   - f_9 * pc_y[k] * fph1_566[k];
    }

#pragma omp simd aligned(t_819, t_820, t_821, t_822, pa_y, pc_y, fph0_567, fph0_568, fph0_569, \
                         fph0_570, fpg_405, fpg_406, fph1_567, fph1_568, fph1_569, \
                         fph1_570 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_819[k] = pa_y[k] * fph0_567[k]
                   - f_9 * pc_y[k] * fph1_567[k];

        t_820[k] = pa_y[k] * fph0_568[k]
                   + f_1 * fpg_405[k]
                   - f_9 * pc_y[k] * fph1_568[k];

        t_821[k] = pa_y[k] * fph0_569[k]
                   - f_9 * pc_y[k] * fph1_569[k];

        t_822[k] = pa_y[k] * fph0_570[k]
                   + f_10 * fpg_406[k]
                   - f_9 * pc_y[k] * fph1_570[k];
    }

#pragma omp simd aligned(t_823, t_824, t_825, pa_y, pc_y, fph0_571, fph0_572, fph0_573, \
                         fpg_407, fpg_408, fph1_571, fph1_572, \
                         fph1_573 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_823[k] = pa_y[k] * fph0_571[k]
                   + f_1 * fpg_407[k]
                   - f_9 * pc_y[k] * fph1_571[k];

        t_824[k] = pa_y[k] * fph0_572[k]
                   - f_9 * pc_y[k] * fph1_572[k];

        t_825[k] = pa_y[k] * fph0_573[k]
                   + f_11 * fpg_408[k]
                   - f_9 * pc_y[k] * fph1_573[k];
    }

#pragma omp simd aligned(t_826, t_827, t_828, pa_y, pc_y, fph0_574, fph0_575, fph0_576, \
                         fpg_409, fpg_410, fph1_574, fph1_575, \
                         fph1_576 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_826[k] = pa_y[k] * fph0_574[k]
                   + f_10 * fpg_409[k]
                   - f_9 * pc_y[k] * fph1_574[k];

        t_827[k] = pa_y[k] * fph0_575[k]
                   + f_1 * fpg_410[k]
                   - f_9 * pc_y[k] * fph1_575[k];

        t_828[k] = pa_y[k] * fph0_576[k]
                   - f_9 * pc_y[k] * fph1_576[k];
    }

#pragma omp simd aligned(t_829, t_830, t_831, t_832, t_833, pc_x, gsg_205, gsg_206, gsg_207, \
                         gsg_208, gsg_209, gpg_595, gpg_596, gpg_597, gpg_598, \
                         gpg_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_829[k] = f_1 * gsg_205[k]
                   + f_4 * pc_x[k] * gpg_595[k];

        t_830[k] = f_1 * gsg_206[k]
                   + f_4 * pc_x[k] * gpg_596[k];

        t_831[k] = f_1 * gsg_207[k]
                   + f_4 * pc_x[k] * gpg_597[k];

        t_832[k] = f_1 * gsg_208[k]
                   + f_4 * pc_x[k] * gpg_598[k];

        t_833[k] = f_1 * gsg_209[k]
                   + f_4 * pc_x[k] * gpg_599[k];
    }
}

static auto
compute_prim_gph_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t fph0, const size_t fpg,
                                                          const size_t fph1, const size_t gsh0,
                                                          const size_t gsg, const size_t gsh1,
                                                          const size_t gpf0, const size_t gpf1,
                                                          const size_t gpg, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = 0.5 / q;
    const auto f_2 = 2.0 / gamma;
    const auto f_3 = 2.0 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = gamma / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 1.5 / q;
    const auto f_12 = 2.5 / q;
    const auto f_15 = 1.5 / gamma;
    const auto f_16 = 1.5 * p / (gamma * q);

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *fph0_587 = buffer.data(fph0 + 587);
    const auto *fph0_609 = buffer.data(fph0 + 609);
    const auto *fph0_611 = buffer.data(fph0 + 611);
    const auto *fph0_614 = buffer.data(fph0 + 614);
    const auto *fph0_618 = buffer.data(fph0 + 618);
    const auto *fph0_624 = buffer.data(fph0 + 624);
    const auto *fph0_626 = buffer.data(fph0 + 626);
    const auto *fph0_627 = buffer.data(fph0 + 627);
    const auto *fph0_629 = buffer.data(fph0 + 629);

    const auto *fpg_370 = buffer.data(fpg + 370);
    const auto *fpg_385 = buffer.data(fpg + 385);
    const auto *fpg_389 = buffer.data(fpg + 389);
    const auto *fpg_400 = buffer.data(fpg + 400);
    const auto *fpg_419 = buffer.data(fpg + 419);
    const auto *fpg_430 = buffer.data(fpg + 430);
    const auto *fpg_432 = buffer.data(fpg + 432);
    const auto *fpg_433 = buffer.data(fpg + 433);
    const auto *fpg_434 = buffer.data(fpg + 434);
    const auto *fpg_445 = buffer.data(fpg + 445);
    const auto *fpg_447 = buffer.data(fpg + 447);
    const auto *fpg_448 = buffer.data(fpg + 448);
    const auto *fpg_449 = buffer.data(fpg + 449);

    const auto *fph1_587 = buffer.data(fph1 + 587);
    const auto *fph1_609 = buffer.data(fph1 + 609);
    const auto *fph1_611 = buffer.data(fph1 + 611);
    const auto *fph1_614 = buffer.data(fph1 + 614);
    const auto *fph1_618 = buffer.data(fph1 + 618);
    const auto *fph1_624 = buffer.data(fph1 + 624);
    const auto *fph1_626 = buffer.data(fph1 + 626);
    const auto *fph1_627 = buffer.data(fph1 + 627);
    const auto *fph1_629 = buffer.data(fph1 + 629);

    const auto *gsh0_288 = buffer.data(gsh0 + 288);
    const auto *gsh0_290 = buffer.data(gsh0 + 290);
    const auto *gsh0_291 = buffer.data(gsh0 + 291);
    const auto *gsh0_294 = buffer.data(gsh0 + 294);
    const auto *gsh0_296 = buffer.data(gsh0 + 296);
    const auto *gsh0_297 = buffer.data(gsh0 + 297);
    const auto *gsh0_299 = buffer.data(gsh0 + 299);
    const auto *gsh0_300 = buffer.data(gsh0 + 300);
    const auto *gsh0_301 = buffer.data(gsh0 + 301);
    const auto *gsh0_303 = buffer.data(gsh0 + 303);
    const auto *gsh0_309 = buffer.data(gsh0 + 309);
    const auto *gsh0_310 = buffer.data(gsh0 + 310);
    const auto *gsh0_311 = buffer.data(gsh0 + 311);
    const auto *gsh0_312 = buffer.data(gsh0 + 312);
    const auto *gsh0_314 = buffer.data(gsh0 + 314);

    const auto *gsg_205 = buffer.data(gsg + 205);
    const auto *gsg_207 = buffer.data(gsg + 207);
    const auto *gsg_208 = buffer.data(gsg + 208);
    const auto *gsg_209 = buffer.data(gsg + 209);
    const auto *gsg_210 = buffer.data(gsg + 210);
    const auto *gsg_212 = buffer.data(gsg + 212);
    const auto *gsg_213 = buffer.data(gsg + 213);
    const auto *gsg_215 = buffer.data(gsg + 215);
    const auto *gsg_216 = buffer.data(gsg + 216);
    const auto *gsg_217 = buffer.data(gsg + 217);
    const auto *gsg_219 = buffer.data(gsg + 219);
    const auto *gsg_220 = buffer.data(gsg + 220);
    const auto *gsg_221 = buffer.data(gsg + 221);
    const auto *gsg_222 = buffer.data(gsg + 222);
    const auto *gsg_223 = buffer.data(gsg + 223);
    const auto *gsg_224 = buffer.data(gsg + 224);

    const auto *gsh1_288 = buffer.data(gsh1 + 288);
    const auto *gsh1_290 = buffer.data(gsh1 + 290);
    const auto *gsh1_291 = buffer.data(gsh1 + 291);
    const auto *gsh1_294 = buffer.data(gsh1 + 294);
    const auto *gsh1_296 = buffer.data(gsh1 + 296);
    const auto *gsh1_297 = buffer.data(gsh1 + 297);
    const auto *gsh1_299 = buffer.data(gsh1 + 299);
    const auto *gsh1_300 = buffer.data(gsh1 + 300);
    const auto *gsh1_301 = buffer.data(gsh1 + 301);
    const auto *gsh1_303 = buffer.data(gsh1 + 303);
    const auto *gsh1_309 = buffer.data(gsh1 + 309);
    const auto *gsh1_310 = buffer.data(gsh1 + 310);
    const auto *gsh1_311 = buffer.data(gsh1 + 311);
    const auto *gsh1_312 = buffer.data(gsh1 + 312);
    const auto *gsh1_314 = buffer.data(gsh1 + 314);

    const auto *gpf0_400 = buffer.data(gpf0 + 400);
    const auto *gpf0_401 = buffer.data(gpf0 + 401);
    const auto *gpf0_402 = buffer.data(gpf0 + 402);
    const auto *gpf0_403 = buffer.data(gpf0 + 403);
    const auto *gpf0_404 = buffer.data(gpf0 + 404);
    const auto *gpf0_405 = buffer.data(gpf0 + 405);
    const auto *gpf0_406 = buffer.data(gpf0 + 406);
    const auto *gpf0_407 = buffer.data(gpf0 + 407);
    const auto *gpf0_408 = buffer.data(gpf0 + 408);
    const auto *gpf0_409 = buffer.data(gpf0 + 409);
    const auto *gpf0_411 = buffer.data(gpf0 + 411);
    const auto *gpf0_413 = buffer.data(gpf0 + 413);
    const auto *gpf0_414 = buffer.data(gpf0 + 414);
    const auto *gpf0_416 = buffer.data(gpf0 + 416);
    const auto *gpf0_417 = buffer.data(gpf0 + 417);
    const auto *gpf0_418 = buffer.data(gpf0 + 418);
    const auto *gpf0_433 = buffer.data(gpf0 + 433);
    const auto *gpf0_436 = buffer.data(gpf0 + 436);
    const auto *gpf0_437 = buffer.data(gpf0 + 437);
    const auto *gpf0_440 = buffer.data(gpf0 + 440);
    const auto *gpf0_442 = buffer.data(gpf0 + 442);
    const auto *gpf0_443 = buffer.data(gpf0 + 443);
    const auto *gpf0_445 = buffer.data(gpf0 + 445);
    const auto *gpf0_446 = buffer.data(gpf0 + 446);
    const auto *gpf0_447 = buffer.data(gpf0 + 447);
    const auto *gpf0_448 = buffer.data(gpf0 + 448);
    const auto *gpf0_449 = buffer.data(gpf0 + 449);

    const auto *gpf1_400 = buffer.data(gpf1 + 400);
    const auto *gpf1_401 = buffer.data(gpf1 + 401);
    const auto *gpf1_402 = buffer.data(gpf1 + 402);
    const auto *gpf1_403 = buffer.data(gpf1 + 403);
    const auto *gpf1_404 = buffer.data(gpf1 + 404);
    const auto *gpf1_405 = buffer.data(gpf1 + 405);
    const auto *gpf1_406 = buffer.data(gpf1 + 406);
    const auto *gpf1_407 = buffer.data(gpf1 + 407);
    const auto *gpf1_408 = buffer.data(gpf1 + 408);
    const auto *gpf1_409 = buffer.data(gpf1 + 409);
    const auto *gpf1_411 = buffer.data(gpf1 + 411);
    const auto *gpf1_413 = buffer.data(gpf1 + 413);
    const auto *gpf1_414 = buffer.data(gpf1 + 414);
    const auto *gpf1_416 = buffer.data(gpf1 + 416);
    const auto *gpf1_417 = buffer.data(gpf1 + 417);
    const auto *gpf1_418 = buffer.data(gpf1 + 418);
    const auto *gpf1_433 = buffer.data(gpf1 + 433);
    const auto *gpf1_436 = buffer.data(gpf1 + 436);
    const auto *gpf1_437 = buffer.data(gpf1 + 437);
    const auto *gpf1_440 = buffer.data(gpf1 + 440);
    const auto *gpf1_442 = buffer.data(gpf1 + 442);
    const auto *gpf1_443 = buffer.data(gpf1 + 443);
    const auto *gpf1_445 = buffer.data(gpf1 + 445);
    const auto *gpf1_446 = buffer.data(gpf1 + 446);
    const auto *gpf1_447 = buffer.data(gpf1 + 447);
    const auto *gpf1_448 = buffer.data(gpf1 + 448);
    const auto *gpf1_449 = buffer.data(gpf1 + 449);

    const auto *gpg_595 = buffer.data(gpg + 595);
    const auto *gpg_599 = buffer.data(gpg + 599);
    const auto *gpg_600 = buffer.data(gpg + 600);
    const auto *gpg_601 = buffer.data(gpg + 601);
    const auto *gpg_602 = buffer.data(gpg + 602);
    const auto *gpg_603 = buffer.data(gpg + 603);
    const auto *gpg_604 = buffer.data(gpg + 604);
    const auto *gpg_605 = buffer.data(gpg + 605);
    const auto *gpg_606 = buffer.data(gpg + 606);
    const auto *gpg_607 = buffer.data(gpg + 607);
    const auto *gpg_608 = buffer.data(gpg + 608);
    const auto *gpg_609 = buffer.data(gpg + 609);
    const auto *gpg_610 = buffer.data(gpg + 610);
    const auto *gpg_611 = buffer.data(gpg + 611);
    const auto *gpg_612 = buffer.data(gpg + 612);
    const auto *gpg_613 = buffer.data(gpg + 613);
    const auto *gpg_614 = buffer.data(gpg + 614);
    const auto *gpg_616 = buffer.data(gpg + 616);
    const auto *gpg_618 = buffer.data(gpg + 618);
    const auto *gpg_619 = buffer.data(gpg + 619);
    const auto *gpg_621 = buffer.data(gpg + 621);
    const auto *gpg_622 = buffer.data(gpg + 622);
    const auto *gpg_623 = buffer.data(gpg + 623);
    const auto *gpg_625 = buffer.data(gpg + 625);
    const auto *gpg_626 = buffer.data(gpg + 626);
    const auto *gpg_627 = buffer.data(gpg + 627);
    const auto *gpg_628 = buffer.data(gpg + 628);
    const auto *gpg_629 = buffer.data(gpg + 629);
    const auto *gpg_630 = buffer.data(gpg + 630);
    const auto *gpg_632 = buffer.data(gpg + 632);
    const auto *gpg_635 = buffer.data(gpg + 635);
    const auto *gpg_640 = buffer.data(gpg + 640);
    const auto *gpg_641 = buffer.data(gpg + 641);
    const auto *gpg_642 = buffer.data(gpg + 642);
    const auto *gpg_643 = buffer.data(gpg + 643);
    const auto *gpg_644 = buffer.data(gpg + 644);
    const auto *gpg_645 = buffer.data(gpg + 645);
    const auto *gpg_647 = buffer.data(gpg + 647);
    const auto *gpg_648 = buffer.data(gpg + 648);
    const auto *gpg_650 = buffer.data(gpg + 650);
    const auto *gpg_651 = buffer.data(gpg + 651);
    const auto *gpg_652 = buffer.data(gpg + 652);
    const auto *gpg_655 = buffer.data(gpg + 655);
    const auto *gpg_656 = buffer.data(gpg + 656);
    const auto *gpg_657 = buffer.data(gpg + 657);
    const auto *gpg_658 = buffer.data(gpg + 658);
    const auto *gpg_659 = buffer.data(gpg + 659);
    const auto *gpg_660 = buffer.data(gpg + 660);
    const auto *gpg_662 = buffer.data(gpg + 662);
    const auto *gpg_663 = buffer.data(gpg + 663);
    const auto *gpg_665 = buffer.data(gpg + 665);
    const auto *gpg_666 = buffer.data(gpg + 666);
    const auto *gpg_667 = buffer.data(gpg + 667);
    const auto *gpg_669 = buffer.data(gpg + 669);
    const auto *gpg_670 = buffer.data(gpg + 670);
    const auto *gpg_671 = buffer.data(gpg + 671);
    const auto *gpg_672 = buffer.data(gpg + 672);
    const auto *gpg_673 = buffer.data(gpg + 673);
    const auto *gpg_674 = buffer.data(gpg + 674);

#pragma omp simd aligned(t_834, t_835, t_836, t_837, pb_x, pc_x, pc_z, fpg_370, gsh0_288, \
                         gsh0_290, gsh0_291, gsh1_288, gsh1_290, gsh1_291, \
                         gpg_595 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_834[k] = pb_x[k] * gsh0_288[k]
                   - f_9 * pc_x[k] * gsh1_288[k];

        t_835[k] = f_11 * fpg_370[k]
                   + f_4 * pc_z[k] * gpg_595[k];

        t_836[k] = pb_x[k] * gsh0_290[k]
                   - f_9 * pc_x[k] * gsh1_290[k];

        t_837[k] = pb_x[k] * gsh0_291[k]
                   - f_9 * pc_x[k] * gsh1_291[k];
    }

#pragma omp simd aligned(t_838, t_839, t_840, pa_y, pc_x, pc_y, fph0_587, fpg_419, fph1_587, \
                         gpf0_400, gpf1_400, gpg_599, gpg_600 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_838[k] = f_1 * fpg_419[k]
                   + f_4 * pc_y[k] * gpg_599[k];

        t_839[k] = pa_y[k] * fph0_587[k]
                   - f_9 * pc_y[k] * fph1_587[k];

        t_840[k] = f_2 * gpf0_400[k]
                   - f_3 * gpf1_400[k]
                   + f_4 * pc_x[k] * gpg_600[k];
    }

#pragma omp simd aligned(t_841, t_842, t_843, pc_x, gpf0_401, gpf0_402, gpf0_403, gpf1_401, \
                         gpf1_402, gpf1_403, gpg_601, gpg_602, \
                         gpg_603 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_841[k] = f_15 * gpf0_401[k]
                   - f_16 * gpf1_401[k]
                   + f_4 * pc_x[k] * gpg_601[k];

        t_842[k] = f_15 * gpf0_402[k]
                   - f_16 * gpf1_402[k]
                   + f_4 * pc_x[k] * gpg_602[k];

        t_843[k] = f_7 * gpf0_403[k]
                   - f_8 * gpf1_403[k]
                   + f_4 * pc_x[k] * gpg_603[k];
    }

#pragma omp simd aligned(t_844, t_845, t_846, pc_x, gpf0_404, gpf0_405, gpf0_406, gpf1_404, \
                         gpf1_405, gpf1_406, gpg_604, gpg_605, \
                         gpg_606 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_844[k] = f_7 * gpf0_404[k]
                   - f_8 * gpf1_404[k]
                   + f_4 * pc_x[k] * gpg_604[k];

        t_845[k] = f_7 * gpf0_405[k]
                   - f_8 * gpf1_405[k]
                   + f_4 * pc_x[k] * gpg_605[k];

        t_846[k] = f_5 * gpf0_406[k]
                   - f_6 * gpf1_406[k]
                   + f_4 * pc_x[k] * gpg_606[k];
    }

#pragma omp simd aligned(t_847, t_848, t_849, t_850, pc_x, gpf0_407, gpf0_408, gpf0_409, \
                         gpf1_407, gpf1_408, gpf1_409, gpg_607, gpg_608, gpg_609, \
                         gpg_610 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_847[k] = f_5 * gpf0_407[k]
                   - f_6 * gpf1_407[k]
                   + f_4 * pc_x[k] * gpg_607[k];

        t_848[k] = f_5 * gpf0_408[k]
                   - f_6 * gpf1_408[k]
                   + f_4 * pc_x[k] * gpg_608[k];

        t_849[k] = f_5 * gpf0_409[k]
                   - f_6 * gpf1_409[k]
                   + f_4 * pc_x[k] * gpg_609[k];

        t_850[k] = f_4 * pc_x[k] * gpg_610[k];
    }

#pragma omp simd aligned(t_851, t_852, t_853, t_854, t_855, pc_x, pc_y, fpg_430, gsg_205, \
                         gpf0_406, gpf1_406, gpg_610, gpg_611, gpg_612, gpg_613, \
                         gpg_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_851[k] = f_4 * pc_x[k] * gpg_611[k];

        t_852[k] = f_4 * pc_x[k] * gpg_612[k];

        t_853[k] = f_4 * pc_x[k] * gpg_613[k];

        t_854[k] = f_4 * pc_x[k] * gpg_614[k];

        t_855[k] = f_1 * fpg_430[k]
                   + f_1 * gsg_205[k]
                   + f_2 * gpf0_406[k]
                   - f_3 * gpf1_406[k]
                   + f_4 * pc_y[k] * gpg_610[k];
    }

#pragma omp simd aligned(t_856, t_857, pc_y, pc_z, fpg_385, fpg_432, gsg_207, gpf0_408, \
                         gpf1_408, gpg_610, gpg_612 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_856[k] = f_11 * fpg_385[k]
                   + f_4 * pc_z[k] * gpg_610[k];

        t_857[k] = f_1 * fpg_432[k]
                   + f_1 * gsg_207[k]
                   + f_7 * gpf0_408[k]
                   - f_8 * gpf1_408[k]
                   + f_4 * pc_y[k] * gpg_612[k];
    }

#pragma omp simd aligned(t_858, t_859, t_860, pc_y, pc_z, fpg_389, fpg_433, fpg_434, gsg_208, \
                         gsg_209, gpf0_409, gpf1_409, gpg_613, \
                         gpg_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_858[k] = f_1 * fpg_433[k]
                   + f_1 * gsg_208[k]
                   + f_5 * gpf0_409[k]
                   - f_6 * gpf1_409[k]
                   + f_4 * pc_y[k] * gpg_613[k];

        t_859[k] = f_1 * fpg_434[k]
                   + f_1 * gsg_209[k]
                   + f_4 * pc_y[k] * gpg_614[k];

        t_860[k] = f_11 * fpg_389[k]
                   + f_2 * gpf0_409[k]
                   - f_3 * gpf1_409[k]
                   + f_4 * pc_z[k] * gpg_614[k];
    }

#pragma omp simd aligned(t_861, t_862, t_863, pa_y, pc_x, pc_y, fph0_609, fph0_611, fph1_609, \
                         fph1_611, gpf0_411, gpf1_411, gpg_616 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_861[k] = pa_y[k] * fph0_609[k]
                   - f_9 * pc_y[k] * fph1_609[k];

        t_862[k] = f_15 * gpf0_411[k]
                   - f_16 * gpf1_411[k]
                   + f_4 * pc_x[k] * gpg_616[k];

        t_863[k] = pa_y[k] * fph0_611[k]
                   - f_9 * pc_y[k] * fph1_611[k];
    }

#pragma omp simd aligned(t_864, t_865, t_866, pa_y, pc_x, pc_y, fph0_614, fph1_614, gpf0_413, \
                         gpf0_414, gpf1_413, gpf1_414, gpg_618, \
                         gpg_619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_864[k] = f_7 * gpf0_413[k]
                   - f_8 * gpf1_413[k]
                   + f_4 * pc_x[k] * gpg_618[k];

        t_865[k] = f_7 * gpf0_414[k]
                   - f_8 * gpf1_414[k]
                   + f_4 * pc_x[k] * gpg_619[k];

        t_866[k] = pa_y[k] * fph0_614[k]
                   - f_9 * pc_y[k] * fph1_614[k];
    }

#pragma omp simd aligned(t_867, t_868, t_869, pc_x, gpf0_416, gpf0_417, gpf0_418, gpf1_416, \
                         gpf1_417, gpf1_418, gpg_621, gpg_622, \
                         gpg_623 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_867[k] = f_5 * gpf0_416[k]
                   - f_6 * gpf1_416[k]
                   + f_4 * pc_x[k] * gpg_621[k];

        t_868[k] = f_5 * gpf0_417[k]
                   - f_6 * gpf1_417[k]
                   + f_4 * pc_x[k] * gpg_622[k];

        t_869[k] = f_5 * gpf0_418[k]
                   - f_6 * gpf1_418[k]
                   + f_4 * pc_x[k] * gpg_623[k];
    }

#pragma omp simd aligned(t_870, t_871, t_872, t_873, t_874, t_875, pa_y, pc_x, pc_y, fph0_618, \
                         fph1_618, gpg_625, gpg_626, gpg_627, gpg_628, \
                         gpg_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_870[k] = pa_y[k] * fph0_618[k]
                   - f_9 * pc_y[k] * fph1_618[k];

        t_871[k] = f_4 * pc_x[k] * gpg_625[k];

        t_872[k] = f_4 * pc_x[k] * gpg_626[k];

        t_873[k] = f_4 * pc_x[k] * gpg_627[k];

        t_874[k] = f_4 * pc_x[k] * gpg_628[k];

        t_875[k] = f_4 * pc_x[k] * gpg_629[k];
    }

#pragma omp simd aligned(t_876, t_877, t_878, pa_y, pc_y, pc_z, fph0_624, fph0_626, fpg_400, \
                         fpg_445, fpg_447, fph1_624, fph1_626, gsg_205, \
                         gpg_625 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_876[k] = pa_y[k] * fph0_624[k]
                   + f_12 * fpg_445[k]
                   - f_9 * pc_y[k] * fph1_624[k];

        t_877[k] = f_11 * fpg_400[k]
                   + f_1 * gsg_205[k]
                   + f_4 * pc_z[k] * gpg_625[k];

        t_878[k] = pa_y[k] * fph0_626[k]
                   + f_11 * fpg_447[k]
                   - f_9 * pc_y[k] * fph1_626[k];
    }

#pragma omp simd aligned(t_879, t_880, t_881, pa_y, pc_y, fph0_627, fph0_629, fpg_448, \
                         fpg_449, fph1_627, fph1_629, gpg_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_879[k] = pa_y[k] * fph0_627[k]
                   + f_10 * fpg_448[k]
                   - f_9 * pc_y[k] * fph1_627[k];

        t_880[k] = f_1 * fpg_449[k]
                   + f_4 * pc_y[k] * gpg_629[k];

        t_881[k] = pa_y[k] * fph0_629[k]
                   - f_9 * pc_y[k] * fph1_629[k];
    }

#pragma omp simd aligned(t_882, t_883, t_884, pb_x, pc_x, pc_y, gsh0_294, gsh0_296, gsg_210, \
                         gsg_212, gsh1_294, gsh1_296, gpg_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_882[k] = pb_x[k] * gsh0_294[k]
                   + f_12 * gsg_210[k]
                   - f_9 * pc_x[k] * gsh1_294[k];

        t_883[k] = f_4 * pc_y[k] * gpg_630[k];

        t_884[k] = pb_x[k] * gsh0_296[k]
                   + f_0 * gsg_212[k]
                   - f_9 * pc_x[k] * gsh1_296[k];
    }

#pragma omp simd aligned(t_885, t_886, t_887, pb_x, pc_x, pc_y, gsh0_297, gsh0_299, gsg_213, \
                         gsg_215, gsh1_297, gsh1_299, gpg_632 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_885[k] = pb_x[k] * gsh0_297[k]
                   + f_11 * gsg_213[k]
                   - f_9 * pc_x[k] * gsh1_297[k];

        t_886[k] = f_4 * pc_y[k] * gpg_632[k];

        t_887[k] = pb_x[k] * gsh0_299[k]
                   + f_11 * gsg_215[k]
                   - f_9 * pc_x[k] * gsh1_299[k];
    }

#pragma omp simd aligned(t_888, t_889, t_890, pb_x, pc_x, pc_y, gsh0_300, gsh0_301, gsg_216, \
                         gsg_217, gsh1_300, gsh1_301, gpg_635 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_888[k] = pb_x[k] * gsh0_300[k]
                   + f_10 * gsg_216[k]
                   - f_9 * pc_x[k] * gsh1_300[k];

        t_889[k] = pb_x[k] * gsh0_301[k]
                   + f_10 * gsg_217[k]
                   - f_9 * pc_x[k] * gsh1_301[k];

        t_890[k] = f_4 * pc_y[k] * gpg_635[k];
    }

#pragma omp simd aligned(t_891, t_892, t_893, t_894, pb_x, pc_x, gsh0_303, gsg_219, gsg_220, \
                         gsg_221, gsg_222, gsh1_303, gpg_640, gpg_641, \
                         gpg_642 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_891[k] = pb_x[k] * gsh0_303[k]
                   + f_10 * gsg_219[k]
                   - f_9 * pc_x[k] * gsh1_303[k];

        t_892[k] = f_1 * gsg_220[k]
                   + f_4 * pc_x[k] * gpg_640[k];

        t_893[k] = f_1 * gsg_221[k]
                   + f_4 * pc_x[k] * gpg_641[k];

        t_894[k] = f_1 * gsg_222[k]
                   + f_4 * pc_x[k] * gpg_642[k];
    }

#pragma omp simd aligned(t_895, t_896, t_897, t_898, pb_x, pc_x, gsh0_309, gsh0_310, gsg_223, \
                         gsg_224, gsh1_309, gsh1_310, gpg_643, \
                         gpg_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_895[k] = f_1 * gsg_223[k]
                   + f_4 * pc_x[k] * gpg_643[k];

        t_896[k] = f_1 * gsg_224[k]
                   + f_4 * pc_x[k] * gpg_644[k];

        t_897[k] = pb_x[k] * gsh0_309[k]
                   - f_9 * pc_x[k] * gsh1_309[k];

        t_898[k] = pb_x[k] * gsh0_310[k]
                   - f_9 * pc_x[k] * gsh1_310[k];
    }

#pragma omp simd aligned(t_899, t_900, t_901, t_902, pb_x, pc_x, pc_y, gsh0_311, gsh0_312, \
                         gsh0_314, gsh1_311, gsh1_312, gsh1_314, \
                         gpg_644 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_899[k] = pb_x[k] * gsh0_311[k]
                   - f_9 * pc_x[k] * gsh1_311[k];

        t_900[k] = pb_x[k] * gsh0_312[k]
                   - f_9 * pc_x[k] * gsh1_312[k];

        t_901[k] = f_4 * pc_y[k] * gpg_644[k];

        t_902[k] = pb_x[k] * gsh0_314[k]
                   - f_9 * pc_x[k] * gsh1_314[k];
    }

#pragma omp simd aligned(t_903, t_904, t_905, t_906, pb_y, pc_x, pc_y, gsh0_294, gsh0_296, \
                         gsg_210, gsh1_294, gsh1_296, gpf0_433, gpf1_433, gpg_645, \
                         gpg_648 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_903[k] = pb_y[k] * gsh0_294[k]
                   - f_9 * pc_y[k] * gsh1_294[k];

        t_904[k] = f_1 * gsg_210[k]
                   + f_4 * pc_y[k] * gpg_645[k];

        t_905[k] = pb_y[k] * gsh0_296[k]
                   - f_9 * pc_y[k] * gsh1_296[k];

        t_906[k] = f_7 * gpf0_433[k]
                   - f_8 * gpf1_433[k]
                   + f_4 * pc_x[k] * gpg_648[k];
    }

#pragma omp simd aligned(t_907, t_908, t_909, pb_y, pc_x, pc_y, gsh0_299, gsg_212, gsh1_299, \
                         gpf0_436, gpf1_436, gpg_647, gpg_651 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_907[k] = f_1 * gsg_212[k]
                   + f_4 * pc_y[k] * gpg_647[k];

        t_908[k] = pb_y[k] * gsh0_299[k]
                   - f_9 * pc_y[k] * gsh1_299[k];

        t_909[k] = f_5 * gpf0_436[k]
                   - f_6 * gpf1_436[k]
                   + f_4 * pc_x[k] * gpg_651[k];
    }

#pragma omp simd aligned(t_910, t_911, t_912, t_913, pb_y, pc_x, pc_y, gsh0_303, gsg_215, \
                         gsh1_303, gpf0_437, gpf1_437, gpg_650, gpg_652, \
                         gpg_655 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_910[k] = f_5 * gpf0_437[k]
                   - f_6 * gpf1_437[k]
                   + f_4 * pc_x[k] * gpg_652[k];

        t_911[k] = f_1 * gsg_215[k]
                   + f_4 * pc_y[k] * gpg_650[k];

        t_912[k] = pb_y[k] * gsh0_303[k]
                   - f_9 * pc_y[k] * gsh1_303[k];

        t_913[k] = f_4 * pc_x[k] * gpg_655[k];
    }

#pragma omp simd aligned(t_914, t_915, t_916, t_917, t_918, pb_y, pc_x, pc_y, gsh0_309, \
                         gsg_220, gsh1_309, gpg_656, gpg_657, gpg_658, \
                         gpg_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_914[k] = f_4 * pc_x[k] * gpg_656[k];

        t_915[k] = f_4 * pc_x[k] * gpg_657[k];

        t_916[k] = f_4 * pc_x[k] * gpg_658[k];

        t_917[k] = f_4 * pc_x[k] * gpg_659[k];

        t_918[k] = pb_y[k] * gsh0_309[k]
                   + f_12 * gsg_220[k]
                   - f_9 * pc_y[k] * gsh1_309[k];
    }

#pragma omp simd aligned(t_919, t_920, t_921, pb_y, pc_y, gsh0_310, gsh0_311, gsh0_312, \
                         gsg_221, gsg_222, gsg_223, gsh1_310, gsh1_311, \
                         gsh1_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_919[k] = pb_y[k] * gsh0_310[k]
                   + f_0 * gsg_221[k]
                   - f_9 * pc_y[k] * gsh1_310[k];

        t_920[k] = pb_y[k] * gsh0_311[k]
                   + f_11 * gsg_222[k]
                   - f_9 * pc_y[k] * gsh1_311[k];

        t_921[k] = pb_y[k] * gsh0_312[k]
                   + f_10 * gsg_223[k]
                   - f_9 * pc_y[k] * gsh1_312[k];
    }

#pragma omp simd aligned(t_922, t_923, t_924, t_925, pb_y, pc_x, pc_y, gsh0_314, gsg_224, \
                         gsh1_314, gpf0_440, gpf1_440, gpg_659, \
                         gpg_660 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_922[k] = f_1 * gsg_224[k]
                   + f_4 * pc_y[k] * gpg_659[k];

        t_923[k] = pb_y[k] * gsh0_314[k]
                   - f_9 * pc_y[k] * gsh1_314[k];

        t_924[k] = f_2 * gpf0_440[k]
                   - f_3 * gpf1_440[k]
                   + f_4 * pc_x[k] * gpg_660[k];

        t_925[k] = f_4 * pc_y[k] * gpg_660[k];
    }

#pragma omp simd aligned(t_926, t_927, t_928, t_929, pc_x, pc_y, gpf0_442, gpf0_443, gpf0_445, \
                         gpf1_442, gpf1_443, gpf1_445, gpg_662, gpg_663, \
                         gpg_665 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_926[k] = f_15 * gpf0_442[k]
                   - f_16 * gpf1_442[k]
                   + f_4 * pc_x[k] * gpg_662[k];

        t_927[k] = f_7 * gpf0_443[k]
                   - f_8 * gpf1_443[k]
                   + f_4 * pc_x[k] * gpg_663[k];

        t_928[k] = f_4 * pc_y[k] * gpg_662[k];

        t_929[k] = f_7 * gpf0_445[k]
                   - f_8 * gpf1_445[k]
                   + f_4 * pc_x[k] * gpg_665[k];
    }

#pragma omp simd aligned(t_930, t_931, t_932, t_933, pc_x, pc_y, gpf0_446, gpf0_447, gpf0_449, \
                         gpf1_446, gpf1_447, gpf1_449, gpg_665, gpg_666, gpg_667, \
                         gpg_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_930[k] = f_5 * gpf0_446[k]
                   - f_6 * gpf1_446[k]
                   + f_4 * pc_x[k] * gpg_666[k];

        t_931[k] = f_5 * gpf0_447[k]
                   - f_6 * gpf1_447[k]
                   + f_4 * pc_x[k] * gpg_667[k];

        t_932[k] = f_4 * pc_y[k] * gpg_665[k];

        t_933[k] = f_5 * gpf0_449[k]
                   - f_6 * gpf1_449[k]
                   + f_4 * pc_x[k] * gpg_669[k];
    }

#pragma omp simd aligned(t_934, t_935, t_936, t_937, t_938, t_939, pc_x, pc_y, gpf0_446, \
                         gpf1_446, gpg_670, gpg_671, gpg_672, gpg_673, \
                         gpg_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_934[k] = f_4 * pc_x[k] * gpg_670[k];

        t_935[k] = f_4 * pc_x[k] * gpg_671[k];

        t_936[k] = f_4 * pc_x[k] * gpg_672[k];

        t_937[k] = f_4 * pc_x[k] * gpg_673[k];

        t_938[k] = f_4 * pc_x[k] * gpg_674[k];

        t_939[k] = f_2 * gpf0_446[k]
                   - f_3 * gpf1_446[k]
                   + f_4 * pc_y[k] * gpg_670[k];
    }

#pragma omp simd aligned(t_940, t_941, t_942, t_943, pc_y, gpf0_447, gpf0_448, gpf0_449, \
                         gpf1_447, gpf1_448, gpf1_449, gpg_671, gpg_672, gpg_673, \
                         gpg_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_940[k] = f_15 * gpf0_447[k]
                   - f_16 * gpf1_447[k]
                   + f_4 * pc_y[k] * gpg_671[k];

        t_941[k] = f_7 * gpf0_448[k]
                   - f_8 * gpf1_448[k]
                   + f_4 * pc_y[k] * gpg_672[k];

        t_942[k] = f_5 * gpf0_449[k]
                   - f_6 * gpf1_449[k]
                   + f_4 * pc_y[k] * gpg_673[k];

        t_943[k] = f_4 * pc_y[k] * gpg_674[k];
    }

#pragma omp simd aligned(t_944, pc_z, fpg_449, gsg_224, gpf0_449, gpf1_449, \
                         gpg_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_944[k] = f_0 * fpg_449[k]
                   + f_1 * gsg_224[k]
                   + f_2 * gpf0_449[k]
                   - f_3 * gpf1_449[k]
                   + f_4 * pc_z[k] * gpg_674[k];
    }
}

auto
compute_prim_gph_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t dph0,
                                                   const size_t dph1, const size_t fph0,
                                                   const size_t fpg, const size_t fph1,
                                                   const size_t gsh0, const size_t gsg,
                                                   const size_t gsh1, const size_t gpf0,
                                                   const size_t gpf1, const size_t gpg,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_gph_three_center_electron_repulsion_0_piece0(buffer, target, pa, pb, pc, dph0,
                                                              dph1, fph0, fpg, fph1, gsh0, gsg,
                                                              gsh1, gpf0, gpf1, gpg, ncols,
                                                              gamma, p, q);

    compute_prim_gph_three_center_electron_repulsion_0_piece1(buffer, target, pa, pb, pc, dph0,
                                                              dph1, fph0, fpg, fph1, gsh0, gsg,
                                                              gsh1, gpf0, gpf1, gpg, ncols,
                                                              gamma, p, q);

    compute_prim_gph_three_center_electron_repulsion_0_piece2(buffer, target, pa, pb, pc, dph0,
                                                              dph1, fph0, fpg, fph1, gsh0, gsg,
                                                              gsh1, gpf0, gpf1, gpg, ncols,
                                                              gamma, p, q);

    compute_prim_gph_three_center_electron_repulsion_0_piece3(buffer, target, pa, pb, pc, dph0,
                                                              dph1, fph0, fpg, fph1, gsh0, gsg,
                                                              gsh1, gpf0, gpf1, gpg, ncols,
                                                              gamma, p, q);

    compute_prim_gph_three_center_electron_repulsion_0_piece4(buffer, target, pa, pb, pc, fph0,
                                                              fpg, fph1, gsh0, gsg, gsh1, gpf0,
                                                              gpf1, gpg, ncols, gamma, p, q);

    compute_prim_gph_three_center_electron_repulsion_0_piece5(buffer, target, pa, pb, pc, fph0,
                                                              fpg, fph1, gsh0, gsg, gsh1, gpf0,
                                                              gpf1, gpg, ncols, gamma, p, q);

    compute_prim_gph_three_center_electron_repulsion_0_piece6(buffer, target, pa, pb, pc, dph0,
                                                              dph1, fph0, fpg, fph1, gsh0, gsg,
                                                              gsh1, gpf0, gpf1, gpg, ncols,
                                                              gamma, p, q);

    compute_prim_gph_three_center_electron_repulsion_0_piece7(buffer, target, pa, pb, pc, fph0,
                                                              fpg, fph1, gsh0, gsg, gsh1, gpf0,
                                                              gpf1, gpg, ncols, gamma, p, q);
}

}  // namespace simdt3ceri
