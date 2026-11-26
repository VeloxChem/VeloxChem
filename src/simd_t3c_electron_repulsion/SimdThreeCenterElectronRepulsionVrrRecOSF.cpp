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


#include "SimdThreeCenterElectronRepulsionVrrRecOSF.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_osf_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsf0,
                                                          const size_t nsd, const size_t nsf1,
                                                          const size_t osp0, const size_t osp1,
                                                          const size_t osd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 5.0 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_9 = 4.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 4.0 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 3.5 / q;
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
    auto *t_124 = buffer.data(target + 124);
    auto *t_125 = buffer.data(target + 125);
    auto *t_126 = buffer.data(target + 126);
    auto *t_127 = buffer.data(target + 127);
    auto *t_128 = buffer.data(target + 128);
    auto *t_129 = buffer.data(target + 129);
    auto *t_130 = buffer.data(target + 130);
    auto *t_131 = buffer.data(target + 131);
    auto *t_132 = buffer.data(target + 132);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsf0_0 = buffer.data(nsf0 + 0);
    const auto *nsf0_6 = buffer.data(nsf0 + 6);
    const auto *nsf0_9 = buffer.data(nsf0 + 9);
    const auto *nsf0_16 = buffer.data(nsf0 + 16);
    const auto *nsf0_20 = buffer.data(nsf0 + 20);
    const auto *nsf0_29 = buffer.data(nsf0 + 29);
    const auto *nsf0_30 = buffer.data(nsf0 + 30);
    const auto *nsf0_36 = buffer.data(nsf0 + 36);
    const auto *nsf0_50 = buffer.data(nsf0 + 50);
    const auto *nsf0_59 = buffer.data(nsf0 + 59);
    const auto *nsf0_60 = buffer.data(nsf0 + 60);
    const auto *nsf0_66 = buffer.data(nsf0 + 66);
    const auto *nsf0_90 = buffer.data(nsf0 + 90);

    const auto *nsd_0 = buffer.data(nsd + 0);
    const auto *nsd_3 = buffer.data(nsd + 3);
    const auto *nsd_5 = buffer.data(nsd + 5);
    const auto *nsd_6 = buffer.data(nsd + 6);
    const auto *nsd_9 = buffer.data(nsd + 9);
    const auto *nsd_11 = buffer.data(nsd + 11);
    const auto *nsd_12 = buffer.data(nsd + 12);
    const auto *nsd_15 = buffer.data(nsd + 15);
    const auto *nsd_17 = buffer.data(nsd + 17);
    const auto *nsd_18 = buffer.data(nsd + 18);
    const auto *nsd_21 = buffer.data(nsd + 21);
    const auto *nsd_23 = buffer.data(nsd + 23);
    const auto *nsd_24 = buffer.data(nsd + 24);
    const auto *nsd_27 = buffer.data(nsd + 27);
    const auto *nsd_28 = buffer.data(nsd + 28);
    const auto *nsd_29 = buffer.data(nsd + 29);
    const auto *nsd_30 = buffer.data(nsd + 30);
    const auto *nsd_33 = buffer.data(nsd + 33);
    const auto *nsd_35 = buffer.data(nsd + 35);
    const auto *nsd_36 = buffer.data(nsd + 36);
    const auto *nsd_39 = buffer.data(nsd + 39);
    const auto *nsd_41 = buffer.data(nsd + 41);
    const auto *nsd_42 = buffer.data(nsd + 42);
    const auto *nsd_45 = buffer.data(nsd + 45);
    const auto *nsd_46 = buffer.data(nsd + 46);
    const auto *nsd_47 = buffer.data(nsd + 47);
    const auto *nsd_48 = buffer.data(nsd + 48);
    const auto *nsd_51 = buffer.data(nsd + 51);
    const auto *nsd_52 = buffer.data(nsd + 52);
    const auto *nsd_53 = buffer.data(nsd + 53);
    const auto *nsd_54 = buffer.data(nsd + 54);
    const auto *nsd_57 = buffer.data(nsd + 57);
    const auto *nsd_59 = buffer.data(nsd + 59);
    const auto *nsd_60 = buffer.data(nsd + 60);
    const auto *nsd_63 = buffer.data(nsd + 63);
    const auto *nsd_65 = buffer.data(nsd + 65);
    const auto *nsd_69 = buffer.data(nsd + 69);
    const auto *nsd_70 = buffer.data(nsd + 70);
    const auto *nsd_71 = buffer.data(nsd + 71);
    const auto *nsd_72 = buffer.data(nsd + 72);
    const auto *nsd_75 = buffer.data(nsd + 75);
    const auto *nsd_76 = buffer.data(nsd + 76);
    const auto *nsd_77 = buffer.data(nsd + 77);

    const auto *nsf1_0 = buffer.data(nsf1 + 0);
    const auto *nsf1_6 = buffer.data(nsf1 + 6);
    const auto *nsf1_9 = buffer.data(nsf1 + 9);
    const auto *nsf1_16 = buffer.data(nsf1 + 16);
    const auto *nsf1_20 = buffer.data(nsf1 + 20);
    const auto *nsf1_29 = buffer.data(nsf1 + 29);
    const auto *nsf1_30 = buffer.data(nsf1 + 30);
    const auto *nsf1_36 = buffer.data(nsf1 + 36);
    const auto *nsf1_50 = buffer.data(nsf1 + 50);
    const auto *nsf1_59 = buffer.data(nsf1 + 59);
    const auto *nsf1_60 = buffer.data(nsf1 + 60);
    const auto *nsf1_66 = buffer.data(nsf1 + 66);
    const auto *nsf1_90 = buffer.data(nsf1 + 90);

    const auto *osp0_0 = buffer.data(osp0 + 0);
    const auto *osp0_1 = buffer.data(osp0 + 1);
    const auto *osp0_2 = buffer.data(osp0 + 2);
    const auto *osp0_4 = buffer.data(osp0 + 4);
    const auto *osp0_8 = buffer.data(osp0 + 8);
    const auto *osp0_9 = buffer.data(osp0 + 9);
    const auto *osp0_10 = buffer.data(osp0 + 10);
    const auto *osp0_11 = buffer.data(osp0 + 11);
    const auto *osp0_15 = buffer.data(osp0 + 15);
    const auto *osp0_16 = buffer.data(osp0 + 16);
    const auto *osp0_17 = buffer.data(osp0 + 17);
    const auto *osp0_18 = buffer.data(osp0 + 18);
    const auto *osp0_19 = buffer.data(osp0 + 19);
    const auto *osp0_20 = buffer.data(osp0 + 20);
    const auto *osp0_23 = buffer.data(osp0 + 23);
    const auto *osp0_25 = buffer.data(osp0 + 25);
    const auto *osp0_27 = buffer.data(osp0 + 27);
    const auto *osp0_28 = buffer.data(osp0 + 28);
    const auto *osp0_29 = buffer.data(osp0 + 29);
    const auto *osp0_30 = buffer.data(osp0 + 30);
    const auto *osp0_31 = buffer.data(osp0 + 31);
    const auto *osp0_32 = buffer.data(osp0 + 32);
    const auto *osp0_35 = buffer.data(osp0 + 35);
    const auto *osp0_36 = buffer.data(osp0 + 36);
    const auto *osp0_37 = buffer.data(osp0 + 37);
    const auto *osp0_38 = buffer.data(osp0 + 38);

    const auto *osp1_0 = buffer.data(osp1 + 0);
    const auto *osp1_1 = buffer.data(osp1 + 1);
    const auto *osp1_2 = buffer.data(osp1 + 2);
    const auto *osp1_4 = buffer.data(osp1 + 4);
    const auto *osp1_8 = buffer.data(osp1 + 8);
    const auto *osp1_9 = buffer.data(osp1 + 9);
    const auto *osp1_10 = buffer.data(osp1 + 10);
    const auto *osp1_11 = buffer.data(osp1 + 11);
    const auto *osp1_15 = buffer.data(osp1 + 15);
    const auto *osp1_16 = buffer.data(osp1 + 16);
    const auto *osp1_17 = buffer.data(osp1 + 17);
    const auto *osp1_18 = buffer.data(osp1 + 18);
    const auto *osp1_19 = buffer.data(osp1 + 19);
    const auto *osp1_20 = buffer.data(osp1 + 20);
    const auto *osp1_23 = buffer.data(osp1 + 23);
    const auto *osp1_25 = buffer.data(osp1 + 25);
    const auto *osp1_27 = buffer.data(osp1 + 27);
    const auto *osp1_28 = buffer.data(osp1 + 28);
    const auto *osp1_29 = buffer.data(osp1 + 29);
    const auto *osp1_30 = buffer.data(osp1 + 30);
    const auto *osp1_31 = buffer.data(osp1 + 31);
    const auto *osp1_32 = buffer.data(osp1 + 32);
    const auto *osp1_35 = buffer.data(osp1 + 35);
    const auto *osp1_36 = buffer.data(osp1 + 36);
    const auto *osp1_37 = buffer.data(osp1 + 37);
    const auto *osp1_38 = buffer.data(osp1 + 38);

    const auto *osd_0 = buffer.data(osd + 0);
    const auto *osd_2 = buffer.data(osd + 2);
    const auto *osd_3 = buffer.data(osd + 3);
    const auto *osd_5 = buffer.data(osd + 5);
    const auto *osd_6 = buffer.data(osd + 6);
    const auto *osd_7 = buffer.data(osd + 7);
    const auto *osd_9 = buffer.data(osd + 9);
    const auto *osd_11 = buffer.data(osd + 11);
    const auto *osd_12 = buffer.data(osd + 12);
    const auto *osd_14 = buffer.data(osd + 14);
    const auto *osd_15 = buffer.data(osd + 15);
    const auto *osd_16 = buffer.data(osd + 16);
    const auto *osd_17 = buffer.data(osd + 17);
    const auto *osd_18 = buffer.data(osd + 18);
    const auto *osd_19 = buffer.data(osd + 19);
    const auto *osd_21 = buffer.data(osd + 21);
    const auto *osd_23 = buffer.data(osd + 23);
    const auto *osd_24 = buffer.data(osd + 24);
    const auto *osd_27 = buffer.data(osd + 27);
    const auto *osd_28 = buffer.data(osd + 28);
    const auto *osd_29 = buffer.data(osd + 29);
    const auto *osd_30 = buffer.data(osd + 30);
    const auto *osd_32 = buffer.data(osd + 32);
    const auto *osd_33 = buffer.data(osd + 33);
    const auto *osd_34 = buffer.data(osd + 34);
    const auto *osd_35 = buffer.data(osd + 35);
    const auto *osd_36 = buffer.data(osd + 36);
    const auto *osd_37 = buffer.data(osd + 37);
    const auto *osd_39 = buffer.data(osd + 39);
    const auto *osd_41 = buffer.data(osd + 41);
    const auto *osd_42 = buffer.data(osd + 42);
    const auto *osd_45 = buffer.data(osd + 45);
    const auto *osd_46 = buffer.data(osd + 46);
    const auto *osd_47 = buffer.data(osd + 47);
    const auto *osd_48 = buffer.data(osd + 48);
    const auto *osd_51 = buffer.data(osd + 51);
    const auto *osd_52 = buffer.data(osd + 52);
    const auto *osd_53 = buffer.data(osd + 53);
    const auto *osd_54 = buffer.data(osd + 54);
    const auto *osd_56 = buffer.data(osd + 56);
    const auto *osd_57 = buffer.data(osd + 57);
    const auto *osd_58 = buffer.data(osd + 58);
    const auto *osd_59 = buffer.data(osd + 59);
    const auto *osd_60 = buffer.data(osd + 60);
    const auto *osd_61 = buffer.data(osd + 61);
    const auto *osd_63 = buffer.data(osd + 63);
    const auto *osd_65 = buffer.data(osd + 65);
    const auto *osd_66 = buffer.data(osd + 66);
    const auto *osd_69 = buffer.data(osd + 69);
    const auto *osd_70 = buffer.data(osd + 70);
    const auto *osd_71 = buffer.data(osd + 71);
    const auto *osd_72 = buffer.data(osd + 72);
    const auto *osd_75 = buffer.data(osd + 75);
    const auto *osd_76 = buffer.data(osd + 76);
    const auto *osd_77 = buffer.data(osd + 77);
    const auto *osd_78 = buffer.data(osd + 78);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, pc_z, nsd_0, nsd_3, osp0_0, \
                         osp1_0, osd_0, osd_2, osd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * nsd_0[k]
                 + f_1 * osp0_0[k]
                 - f_2 * osp1_0[k]
                 + f_3 * pc_x[k] * osd_0[k];

        t_1[k] = f_3 * pc_y[k] * osd_0[k];

        t_2[k] = f_3 * pc_z[k] * osd_0[k];

        t_3[k] = f_0 * nsd_3[k]
                 + f_3 * pc_x[k] * osd_3[k];

        t_4[k] = f_3 * pc_y[k] * osd_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pc_x, pc_y, pc_z, nsd_5, osp0_1, osp0_2, \
                         osp1_1, osp1_2, osd_3, osd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * nsd_5[k]
                 + f_3 * pc_x[k] * osd_5[k];

        t_6[k] = f_1 * osp0_1[k]
                 - f_2 * osp1_1[k]
                 + f_3 * pc_y[k] * osd_3[k];

        t_7[k] = f_3 * pc_z[k] * osd_3[k];

        t_8[k] = f_3 * pc_y[k] * osd_5[k];

        t_9[k] = f_1 * osp0_2[k]
                 - f_2 * osp1_2[k]
                 + f_3 * pc_z[k] * osd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pa_y, pc_x, pc_y, pc_z, nsf0_0, nsd_0, \
                         nsd_9, nsf1_0, osd_6, osd_7, osd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * nsf0_0[k]
                  - f_4 * pc_y[k] * nsf1_0[k];

        t_11[k] = f_5 * nsd_0[k]
                  + f_3 * pc_y[k] * osd_6[k];

        t_12[k] = f_3 * pc_z[k] * osd_6[k];

        t_13[k] = f_6 * nsd_9[k]
                  + f_3 * pc_x[k] * osd_9[k];

        t_14[k] = f_3 * pc_z[k] * osd_7[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pc_x, pc_y, pc_z, nsd_3, nsd_5, nsd_11, \
                         osp0_4, osp1_4, osd_9, osd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_6 * nsd_11[k]
                  + f_3 * pc_x[k] * osd_11[k];

        t_16[k] = f_5 * nsd_3[k]
                  + f_1 * osp0_4[k]
                  - f_2 * osp1_4[k]
                  + f_3 * pc_y[k] * osd_9[k];

        t_17[k] = f_3 * pc_z[k] * osd_9[k];

        t_18[k] = f_5 * nsd_5[k]
                  + f_3 * pc_y[k] * osd_11[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_y, pa_z, pc_y, pc_z, nsf0_0, nsf0_9, \
                         nsd_0, nsf1_0, nsf1_9, osd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_y[k] * nsf0_9[k]
                  - f_4 * pc_y[k] * nsf1_9[k];

        t_20[k] = pa_z[k] * nsf0_0[k]
                  - f_4 * pc_z[k] * nsf1_0[k];

        t_21[k] = f_3 * pc_y[k] * osd_12[k];

        t_22[k] = f_5 * nsd_0[k]
                  + f_3 * pc_z[k] * osd_12[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_z, pc_x, pc_y, pc_z, nsf0_6, nsd_15, \
                         nsd_17, nsf1_6, osd_14, osd_15, osd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_6 * nsd_15[k]
                  + f_3 * pc_x[k] * osd_15[k];

        t_24[k] = f_3 * pc_y[k] * osd_14[k];

        t_25[k] = f_6 * nsd_17[k]
                  + f_3 * pc_x[k] * osd_17[k];

        t_26[k] = pa_z[k] * nsf0_6[k]
                  - f_4 * pc_z[k] * nsf1_6[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pc_x, pc_y, pc_z, nsd_5, nsd_18, osp0_8, \
                         osp0_9, osp1_8, osp1_9, osd_16, osd_17, \
                         osd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_7 * osp0_8[k]
                  - f_8 * osp1_8[k]
                  + f_3 * pc_y[k] * osd_16[k];

        t_28[k] = f_3 * pc_y[k] * osd_17[k];

        t_29[k] = f_5 * nsd_5[k]
                  + f_1 * osp0_8[k]
                  - f_2 * osp1_8[k]
                  + f_3 * pc_z[k] * osd_17[k];

        t_30[k] = f_9 * nsd_18[k]
                  + f_1 * osp0_9[k]
                  - f_2 * osp1_9[k]
                  + f_3 * pc_x[k] * osd_18[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pc_x, pc_y, pc_z, nsd_6, nsd_21, \
                         nsd_23, osd_18, osd_19, osd_21, osd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_10 * nsd_6[k]
                  + f_3 * pc_y[k] * osd_18[k];

        t_32[k] = f_3 * pc_z[k] * osd_18[k];

        t_33[k] = f_9 * nsd_21[k]
                  + f_3 * pc_x[k] * osd_21[k];

        t_34[k] = f_3 * pc_z[k] * osd_19[k];

        t_35[k] = f_9 * nsd_23[k]
                  + f_3 * pc_x[k] * osd_23[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pc_y, pc_z, nsd_9, nsd_11, osp0_10, osp0_11, \
                         osp1_10, osp1_11, osd_21, osd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_10 * nsd_9[k]
                  + f_1 * osp0_10[k]
                  - f_2 * osp1_10[k]
                  + f_3 * pc_y[k] * osd_21[k];

        t_37[k] = f_3 * pc_z[k] * osd_21[k];

        t_38[k] = f_10 * nsd_11[k]
                  + f_3 * pc_y[k] * osd_23[k];

        t_39[k] = f_1 * osp0_11[k]
                  - f_2 * osp1_11[k]
                  + f_3 * pc_z[k] * osd_23[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_y, pc_x, pc_y, pc_z, nsf0_20, nsd_6, \
                         nsd_12, nsd_27, nsf1_20, osd_24, osd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_y[k] * nsf0_20[k]
                  - f_4 * pc_y[k] * nsf1_20[k];

        t_41[k] = f_5 * nsd_12[k]
                  + f_3 * pc_y[k] * osd_24[k];

        t_42[k] = f_5 * nsd_6[k]
                  + f_3 * pc_z[k] * osd_24[k];

        t_43[k] = f_9 * nsd_27[k]
                  + f_3 * pc_x[k] * osd_27[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_z, pc_x, pc_z, nsf0_16, nsd_9, nsd_28, \
                         nsd_29, nsf1_16, osd_27, osd_28, osd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_9 * nsd_28[k]
                  + f_3 * pc_x[k] * osd_28[k];

        t_45[k] = f_9 * nsd_29[k]
                  + f_3 * pc_x[k] * osd_29[k];

        t_46[k] = pa_z[k] * nsf0_16[k]
                  - f_4 * pc_z[k] * nsf1_16[k];

        t_47[k] = f_5 * nsd_9[k]
                  + f_3 * pc_z[k] * osd_27[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_y, pc_x, pc_y, nsf0_29, nsd_17, nsd_30, \
                         nsf1_29, osp0_15, osp1_15, osd_29, osd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_5 * nsd_17[k]
                  + f_3 * pc_y[k] * osd_29[k];

        t_49[k] = pa_y[k] * nsf0_29[k]
                  - f_4 * pc_y[k] * nsf1_29[k];

        t_50[k] = f_9 * nsd_30[k]
                  + f_1 * osp0_15[k]
                  - f_2 * osp1_15[k]
                  + f_3 * pc_x[k] * osd_30[k];

        t_51[k] = f_3 * pc_y[k] * osd_30[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pc_x, pc_y, pc_z, nsd_12, nsd_33, nsd_35, \
                         osd_30, osd_32, osd_33, osd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_10 * nsd_12[k]
                  + f_3 * pc_z[k] * osd_30[k];

        t_53[k] = f_9 * nsd_33[k]
                  + f_3 * pc_x[k] * osd_33[k];

        t_54[k] = f_3 * pc_y[k] * osd_32[k];

        t_55[k] = f_9 * nsd_35[k]
                  + f_3 * pc_x[k] * osd_35[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pc_y, pc_z, nsd_17, osp0_16, osp0_17, \
                         osp1_16, osp1_17, osd_33, osd_34, osd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_1 * osp0_16[k]
                  - f_2 * osp1_16[k]
                  + f_3 * pc_y[k] * osd_33[k];

        t_57[k] = f_7 * osp0_17[k]
                  - f_8 * osp1_17[k]
                  + f_3 * pc_y[k] * osd_34[k];

        t_58[k] = f_3 * pc_y[k] * osd_35[k];

        t_59[k] = f_10 * nsd_17[k]
                  + f_1 * osp0_17[k]
                  - f_2 * osp1_17[k]
                  + f_3 * pc_z[k] * osd_35[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, pc_x, pc_y, pc_z, nsd_18, nsd_36, \
                         nsd_39, osp0_18, osp1_18, osd_36, osd_37, \
                         osd_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_11 * nsd_36[k]
                  + f_1 * osp0_18[k]
                  - f_2 * osp1_18[k]
                  + f_3 * pc_x[k] * osd_36[k];

        t_61[k] = f_12 * nsd_18[k]
                  + f_3 * pc_y[k] * osd_36[k];

        t_62[k] = f_3 * pc_z[k] * osd_36[k];

        t_63[k] = f_11 * nsd_39[k]
                  + f_3 * pc_x[k] * osd_39[k];

        t_64[k] = f_3 * pc_z[k] * osd_37[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pc_x, pc_y, pc_z, nsd_21, nsd_23, nsd_41, \
                         osp0_19, osp1_19, osd_39, osd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_11 * nsd_41[k]
                  + f_3 * pc_x[k] * osd_41[k];

        t_66[k] = f_12 * nsd_21[k]
                  + f_1 * osp0_19[k]
                  - f_2 * osp1_19[k]
                  + f_3 * pc_y[k] * osd_39[k];

        t_67[k] = f_3 * pc_z[k] * osd_39[k];

        t_68[k] = f_12 * nsd_23[k]
                  + f_3 * pc_y[k] * osd_41[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_z, pc_y, pc_z, nsf0_30, nsd_18, nsd_24, \
                         nsf1_30, osp0_20, osp1_20, osd_41, osd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_1 * osp0_20[k]
                  - f_2 * osp1_20[k]
                  + f_3 * pc_z[k] * osd_41[k];

        t_70[k] = pa_z[k] * nsf0_30[k]
                  - f_4 * pc_z[k] * nsf1_30[k];

        t_71[k] = f_10 * nsd_24[k]
                  + f_3 * pc_y[k] * osd_42[k];

        t_72[k] = f_5 * nsd_18[k]
                  + f_3 * pc_z[k] * osd_42[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_z, pc_x, pc_z, nsf0_36, nsd_45, nsd_46, \
                         nsd_47, nsf1_36, osd_45, osd_46, osd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_11 * nsd_45[k]
                  + f_3 * pc_x[k] * osd_45[k];

        t_74[k] = f_11 * nsd_46[k]
                  + f_3 * pc_x[k] * osd_46[k];

        t_75[k] = f_11 * nsd_47[k]
                  + f_3 * pc_x[k] * osd_47[k];

        t_76[k] = pa_z[k] * nsf0_36[k]
                  - f_4 * pc_z[k] * nsf1_36[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_y, pc_y, pc_z, nsf0_50, nsd_21, nsd_23, \
                         nsd_29, nsf1_50, osp0_23, osp1_23, osd_45, \
                         osd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_5 * nsd_21[k]
                  + f_3 * pc_z[k] * osd_45[k];

        t_78[k] = f_10 * nsd_29[k]
                  + f_3 * pc_y[k] * osd_47[k];

        t_79[k] = f_5 * nsd_23[k]
                  + f_1 * osp0_23[k]
                  - f_2 * osp1_23[k]
                  + f_3 * pc_z[k] * osd_47[k];

        t_80[k] = pa_y[k] * nsf0_50[k]
                  - f_4 * pc_y[k] * nsf1_50[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pc_x, pc_y, pc_z, nsd_24, nsd_30, nsd_51, \
                         nsd_52, osd_48, osd_51, osd_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_5 * nsd_30[k]
                  + f_3 * pc_y[k] * osd_48[k];

        t_82[k] = f_10 * nsd_24[k]
                  + f_3 * pc_z[k] * osd_48[k];

        t_83[k] = f_11 * nsd_51[k]
                  + f_3 * pc_x[k] * osd_51[k];

        t_84[k] = f_11 * nsd_52[k]
                  + f_3 * pc_x[k] * osd_52[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pc_x, pc_y, pc_z, nsd_27, nsd_33, nsd_35, \
                         nsd_53, osp0_25, osp1_25, osd_51, osd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = f_11 * nsd_53[k]
                  + f_3 * pc_x[k] * osd_53[k];

        t_86[k] = f_5 * nsd_33[k]
                  + f_1 * osp0_25[k]
                  - f_2 * osp1_25[k]
                  + f_3 * pc_y[k] * osd_51[k];

        t_87[k] = f_10 * nsd_27[k]
                  + f_3 * pc_z[k] * osd_51[k];

        t_88[k] = f_5 * nsd_35[k]
                  + f_3 * pc_y[k] * osd_53[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pa_y, pc_x, pc_y, pc_z, nsf0_59, nsd_30, \
                         nsd_54, nsf1_59, osp0_27, osp1_27, osd_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = pa_y[k] * nsf0_59[k]
                  - f_4 * pc_y[k] * nsf1_59[k];

        t_90[k] = f_11 * nsd_54[k]
                  + f_1 * osp0_27[k]
                  - f_2 * osp1_27[k]
                  + f_3 * pc_x[k] * osd_54[k];

        t_91[k] = f_3 * pc_y[k] * osd_54[k];

        t_92[k] = f_12 * nsd_30[k]
                  + f_3 * pc_z[k] * osd_54[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pc_x, pc_y, nsd_57, nsd_59, osp0_28, osp1_28, \
                         osd_56, osd_57, osd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_11 * nsd_57[k]
                  + f_3 * pc_x[k] * osd_57[k];

        t_94[k] = f_3 * pc_y[k] * osd_56[k];

        t_95[k] = f_11 * nsd_59[k]
                  + f_3 * pc_x[k] * osd_59[k];

        t_96[k] = f_1 * osp0_28[k]
                  - f_2 * osp1_28[k]
                  + f_3 * pc_y[k] * osd_57[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pc_x, pc_y, pc_z, nsd_35, nsd_60, osp0_29, \
                         osp0_30, osp1_29, osp1_30, osd_58, osd_59, \
                         osd_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_7 * osp0_29[k]
                  - f_8 * osp1_29[k]
                  + f_3 * pc_y[k] * osd_58[k];

        t_98[k] = f_3 * pc_y[k] * osd_59[k];

        t_99[k] = f_12 * nsd_35[k]
                  + f_1 * osp0_29[k]
                  - f_2 * osp1_29[k]
                  + f_3 * pc_z[k] * osd_59[k];

        t_100[k] = f_13 * nsd_60[k]
                   + f_1 * osp0_30[k]
                   - f_2 * osp1_30[k]
                   + f_3 * pc_x[k] * osd_60[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, pc_x, pc_y, pc_z, nsd_36, nsd_63, \
                         nsd_65, osd_60, osd_61, osd_63, osd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_14 * nsd_36[k]
                   + f_3 * pc_y[k] * osd_60[k];

        t_102[k] = f_3 * pc_z[k] * osd_60[k];

        t_103[k] = f_13 * nsd_63[k]
                   + f_3 * pc_x[k] * osd_63[k];

        t_104[k] = f_3 * pc_z[k] * osd_61[k];

        t_105[k] = f_13 * nsd_65[k]
                   + f_3 * pc_x[k] * osd_65[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pc_y, pc_z, nsd_39, nsd_41, osp0_31, \
                         osp0_32, osp1_31, osp1_32, osd_63, osd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_14 * nsd_39[k]
                   + f_1 * osp0_31[k]
                   - f_2 * osp1_31[k]
                   + f_3 * pc_y[k] * osd_63[k];

        t_107[k] = f_3 * pc_z[k] * osd_63[k];

        t_108[k] = f_14 * nsd_41[k]
                   + f_3 * pc_y[k] * osd_65[k];

        t_109[k] = f_1 * osp0_32[k]
                   - f_2 * osp1_32[k]
                   + f_3 * pc_z[k] * osd_65[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pa_z, pc_x, pc_y, pc_z, nsf0_60, nsd_36, \
                         nsd_42, nsd_69, nsf1_60, osd_66, osd_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pa_z[k] * nsf0_60[k]
                   - f_4 * pc_z[k] * nsf1_60[k];

        t_111[k] = f_12 * nsd_42[k]
                   + f_3 * pc_y[k] * osd_66[k];

        t_112[k] = f_5 * nsd_36[k]
                   + f_3 * pc_z[k] * osd_66[k];

        t_113[k] = f_13 * nsd_69[k]
                   + f_3 * pc_x[k] * osd_69[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pa_z, pc_x, pc_z, nsf0_66, nsd_39, \
                         nsd_70, nsd_71, nsf1_66, osd_69, osd_70, \
                         osd_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_13 * nsd_70[k]
                   + f_3 * pc_x[k] * osd_70[k];

        t_115[k] = f_13 * nsd_71[k]
                   + f_3 * pc_x[k] * osd_71[k];

        t_116[k] = pa_z[k] * nsf0_66[k]
                   - f_4 * pc_z[k] * nsf1_66[k];

        t_117[k] = f_5 * nsd_39[k]
                   + f_3 * pc_z[k] * osd_69[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pc_x, pc_y, pc_z, nsd_41, nsd_47, nsd_72, \
                         osp0_35, osp0_36, osp1_35, osp1_36, osd_71, \
                         osd_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_12 * nsd_47[k]
                   + f_3 * pc_y[k] * osd_71[k];

        t_119[k] = f_5 * nsd_41[k]
                   + f_1 * osp0_35[k]
                   - f_2 * osp1_35[k]
                   + f_3 * pc_z[k] * osd_71[k];

        t_120[k] = f_13 * nsd_72[k]
                   + f_1 * osp0_36[k]
                   - f_2 * osp1_36[k]
                   + f_3 * pc_x[k] * osd_72[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pc_x, pc_y, pc_z, nsd_42, nsd_48, nsd_75, \
                         nsd_76, osd_72, osd_75, osd_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_10 * nsd_48[k]
                   + f_3 * pc_y[k] * osd_72[k];

        t_122[k] = f_10 * nsd_42[k]
                   + f_3 * pc_z[k] * osd_72[k];

        t_123[k] = f_13 * nsd_75[k]
                   + f_3 * pc_x[k] * osd_75[k];

        t_124[k] = f_13 * nsd_76[k]
                   + f_3 * pc_x[k] * osd_76[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pc_x, pc_y, pc_z, nsd_45, nsd_51, nsd_53, \
                         nsd_77, osp0_37, osp1_37, osd_75, osd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_13 * nsd_77[k]
                   + f_3 * pc_x[k] * osd_77[k];

        t_126[k] = f_10 * nsd_51[k]
                   + f_1 * osp0_37[k]
                   - f_2 * osp1_37[k]
                   + f_3 * pc_y[k] * osd_75[k];

        t_127[k] = f_10 * nsd_45[k]
                   + f_3 * pc_z[k] * osd_75[k];

        t_128[k] = f_10 * nsd_53[k]
                   + f_3 * pc_y[k] * osd_77[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, pa_y, pc_y, pc_z, nsf0_90, nsd_47, \
                         nsd_48, nsd_54, nsf1_90, osp0_38, osp1_38, osd_77, \
                         osd_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_10 * nsd_47[k]
                   + f_1 * osp0_38[k]
                   - f_2 * osp1_38[k]
                   + f_3 * pc_z[k] * osd_77[k];

        t_130[k] = pa_y[k] * nsf0_90[k]
                   - f_4 * pc_y[k] * nsf1_90[k];

        t_131[k] = f_5 * nsd_54[k]
                   + f_3 * pc_y[k] * osd_78[k];

        t_132[k] = f_12 * nsd_48[k]
                   + f_3 * pc_z[k] * osd_78[k];
    }
}

static auto
compute_prim_osf_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsf0,
                                                          const size_t nsd, const size_t nsf1,
                                                          const size_t osp0, const size_t osp1,
                                                          const size_t osd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_10 = 1.0 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 3.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.5 / q;

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
    auto *t_249 = buffer.data(target + 249);
    auto *t_250 = buffer.data(target + 250);
    auto *t_251 = buffer.data(target + 251);
    auto *t_252 = buffer.data(target + 252);
    auto *t_253 = buffer.data(target + 253);
    auto *t_254 = buffer.data(target + 254);
    auto *t_255 = buffer.data(target + 255);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsf0_99 = buffer.data(nsf0 + 99);
    const auto *nsf0_100 = buffer.data(nsf0 + 100);
    const auto *nsf0_106 = buffer.data(nsf0 + 106);
    const auto *nsf0_140 = buffer.data(nsf0 + 140);
    const auto *nsf0_149 = buffer.data(nsf0 + 149);
    const auto *nsf0_150 = buffer.data(nsf0 + 150);
    const auto *nsf0_156 = buffer.data(nsf0 + 156);

    const auto *nsd_51 = buffer.data(nsd + 51);
    const auto *nsd_54 = buffer.data(nsd + 54);
    const auto *nsd_57 = buffer.data(nsd + 57);
    const auto *nsd_59 = buffer.data(nsd + 59);
    const auto *nsd_60 = buffer.data(nsd + 60);
    const auto *nsd_63 = buffer.data(nsd + 63);
    const auto *nsd_65 = buffer.data(nsd + 65);
    const auto *nsd_66 = buffer.data(nsd + 66);
    const auto *nsd_69 = buffer.data(nsd + 69);
    const auto *nsd_71 = buffer.data(nsd + 71);
    const auto *nsd_72 = buffer.data(nsd + 72);
    const auto *nsd_75 = buffer.data(nsd + 75);
    const auto *nsd_77 = buffer.data(nsd + 77);
    const auto *nsd_78 = buffer.data(nsd + 78);
    const auto *nsd_81 = buffer.data(nsd + 81);
    const auto *nsd_82 = buffer.data(nsd + 82);
    const auto *nsd_83 = buffer.data(nsd + 83);
    const auto *nsd_84 = buffer.data(nsd + 84);
    const auto *nsd_87 = buffer.data(nsd + 87);
    const auto *nsd_89 = buffer.data(nsd + 89);
    const auto *nsd_90 = buffer.data(nsd + 90);
    const auto *nsd_93 = buffer.data(nsd + 93);
    const auto *nsd_95 = buffer.data(nsd + 95);
    const auto *nsd_96 = buffer.data(nsd + 96);
    const auto *nsd_99 = buffer.data(nsd + 99);
    const auto *nsd_100 = buffer.data(nsd + 100);
    const auto *nsd_101 = buffer.data(nsd + 101);
    const auto *nsd_102 = buffer.data(nsd + 102);
    const auto *nsd_105 = buffer.data(nsd + 105);
    const auto *nsd_106 = buffer.data(nsd + 106);
    const auto *nsd_107 = buffer.data(nsd + 107);
    const auto *nsd_108 = buffer.data(nsd + 108);
    const auto *nsd_111 = buffer.data(nsd + 111);
    const auto *nsd_112 = buffer.data(nsd + 112);
    const auto *nsd_113 = buffer.data(nsd + 113);
    const auto *nsd_114 = buffer.data(nsd + 114);
    const auto *nsd_117 = buffer.data(nsd + 117);
    const auto *nsd_118 = buffer.data(nsd + 118);
    const auto *nsd_119 = buffer.data(nsd + 119);
    const auto *nsd_120 = buffer.data(nsd + 120);
    const auto *nsd_123 = buffer.data(nsd + 123);
    const auto *nsd_125 = buffer.data(nsd + 125);
    const auto *nsd_126 = buffer.data(nsd + 126);
    const auto *nsd_129 = buffer.data(nsd + 129);
    const auto *nsd_131 = buffer.data(nsd + 131);
    const auto *nsd_135 = buffer.data(nsd + 135);
    const auto *nsd_136 = buffer.data(nsd + 136);
    const auto *nsd_137 = buffer.data(nsd + 137);
    const auto *nsd_138 = buffer.data(nsd + 138);
    const auto *nsd_141 = buffer.data(nsd + 141);
    const auto *nsd_142 = buffer.data(nsd + 142);
    const auto *nsd_143 = buffer.data(nsd + 143);
    const auto *nsd_144 = buffer.data(nsd + 144);
    const auto *nsd_147 = buffer.data(nsd + 147);
    const auto *nsd_148 = buffer.data(nsd + 148);
    const auto *nsd_149 = buffer.data(nsd + 149);
    const auto *nsd_150 = buffer.data(nsd + 150);
    const auto *nsd_153 = buffer.data(nsd + 153);
    const auto *nsd_154 = buffer.data(nsd + 154);
    const auto *nsd_155 = buffer.data(nsd + 155);

    const auto *nsf1_99 = buffer.data(nsf1 + 99);
    const auto *nsf1_100 = buffer.data(nsf1 + 100);
    const auto *nsf1_106 = buffer.data(nsf1 + 106);
    const auto *nsf1_140 = buffer.data(nsf1 + 140);
    const auto *nsf1_149 = buffer.data(nsf1 + 149);
    const auto *nsf1_150 = buffer.data(nsf1 + 150);
    const auto *nsf1_156 = buffer.data(nsf1 + 156);

    const auto *osp0_40 = buffer.data(osp0 + 40);
    const auto *osp0_42 = buffer.data(osp0 + 42);
    const auto *osp0_43 = buffer.data(osp0 + 43);
    const auto *osp0_44 = buffer.data(osp0 + 44);
    const auto *osp0_45 = buffer.data(osp0 + 45);
    const auto *osp0_46 = buffer.data(osp0 + 46);
    const auto *osp0_47 = buffer.data(osp0 + 47);
    const auto *osp0_50 = buffer.data(osp0 + 50);
    const auto *osp0_51 = buffer.data(osp0 + 51);
    const auto *osp0_52 = buffer.data(osp0 + 52);
    const auto *osp0_53 = buffer.data(osp0 + 53);
    const auto *osp0_54 = buffer.data(osp0 + 54);
    const auto *osp0_55 = buffer.data(osp0 + 55);
    const auto *osp0_56 = buffer.data(osp0 + 56);
    const auto *osp0_58 = buffer.data(osp0 + 58);
    const auto *osp0_60 = buffer.data(osp0 + 60);
    const auto *osp0_61 = buffer.data(osp0 + 61);
    const auto *osp0_62 = buffer.data(osp0 + 62);
    const auto *osp0_63 = buffer.data(osp0 + 63);
    const auto *osp0_64 = buffer.data(osp0 + 64);
    const auto *osp0_65 = buffer.data(osp0 + 65);
    const auto *osp0_68 = buffer.data(osp0 + 68);
    const auto *osp0_69 = buffer.data(osp0 + 69);
    const auto *osp0_70 = buffer.data(osp0 + 70);
    const auto *osp0_71 = buffer.data(osp0 + 71);
    const auto *osp0_72 = buffer.data(osp0 + 72);
    const auto *osp0_73 = buffer.data(osp0 + 73);
    const auto *osp0_74 = buffer.data(osp0 + 74);
    const auto *osp0_75 = buffer.data(osp0 + 75);

    const auto *osp1_40 = buffer.data(osp1 + 40);
    const auto *osp1_42 = buffer.data(osp1 + 42);
    const auto *osp1_43 = buffer.data(osp1 + 43);
    const auto *osp1_44 = buffer.data(osp1 + 44);
    const auto *osp1_45 = buffer.data(osp1 + 45);
    const auto *osp1_46 = buffer.data(osp1 + 46);
    const auto *osp1_47 = buffer.data(osp1 + 47);
    const auto *osp1_50 = buffer.data(osp1 + 50);
    const auto *osp1_51 = buffer.data(osp1 + 51);
    const auto *osp1_52 = buffer.data(osp1 + 52);
    const auto *osp1_53 = buffer.data(osp1 + 53);
    const auto *osp1_54 = buffer.data(osp1 + 54);
    const auto *osp1_55 = buffer.data(osp1 + 55);
    const auto *osp1_56 = buffer.data(osp1 + 56);
    const auto *osp1_58 = buffer.data(osp1 + 58);
    const auto *osp1_60 = buffer.data(osp1 + 60);
    const auto *osp1_61 = buffer.data(osp1 + 61);
    const auto *osp1_62 = buffer.data(osp1 + 62);
    const auto *osp1_63 = buffer.data(osp1 + 63);
    const auto *osp1_64 = buffer.data(osp1 + 64);
    const auto *osp1_65 = buffer.data(osp1 + 65);
    const auto *osp1_68 = buffer.data(osp1 + 68);
    const auto *osp1_69 = buffer.data(osp1 + 69);
    const auto *osp1_70 = buffer.data(osp1 + 70);
    const auto *osp1_71 = buffer.data(osp1 + 71);
    const auto *osp1_72 = buffer.data(osp1 + 72);
    const auto *osp1_73 = buffer.data(osp1 + 73);
    const auto *osp1_74 = buffer.data(osp1 + 74);
    const auto *osp1_75 = buffer.data(osp1 + 75);

    const auto *osd_81 = buffer.data(osd + 81);
    const auto *osd_82 = buffer.data(osd + 82);
    const auto *osd_83 = buffer.data(osd + 83);
    const auto *osd_84 = buffer.data(osd + 84);
    const auto *osd_86 = buffer.data(osd + 86);
    const auto *osd_87 = buffer.data(osd + 87);
    const auto *osd_88 = buffer.data(osd + 88);
    const auto *osd_89 = buffer.data(osd + 89);
    const auto *osd_90 = buffer.data(osd + 90);
    const auto *osd_91 = buffer.data(osd + 91);
    const auto *osd_93 = buffer.data(osd + 93);
    const auto *osd_95 = buffer.data(osd + 95);
    const auto *osd_96 = buffer.data(osd + 96);
    const auto *osd_99 = buffer.data(osd + 99);
    const auto *osd_100 = buffer.data(osd + 100);
    const auto *osd_101 = buffer.data(osd + 101);
    const auto *osd_102 = buffer.data(osd + 102);
    const auto *osd_105 = buffer.data(osd + 105);
    const auto *osd_106 = buffer.data(osd + 106);
    const auto *osd_107 = buffer.data(osd + 107);
    const auto *osd_108 = buffer.data(osd + 108);
    const auto *osd_111 = buffer.data(osd + 111);
    const auto *osd_112 = buffer.data(osd + 112);
    const auto *osd_113 = buffer.data(osd + 113);
    const auto *osd_114 = buffer.data(osd + 114);
    const auto *osd_117 = buffer.data(osd + 117);
    const auto *osd_118 = buffer.data(osd + 118);
    const auto *osd_119 = buffer.data(osd + 119);
    const auto *osd_120 = buffer.data(osd + 120);
    const auto *osd_122 = buffer.data(osd + 122);
    const auto *osd_123 = buffer.data(osd + 123);
    const auto *osd_124 = buffer.data(osd + 124);
    const auto *osd_125 = buffer.data(osd + 125);
    const auto *osd_126 = buffer.data(osd + 126);
    const auto *osd_127 = buffer.data(osd + 127);
    const auto *osd_129 = buffer.data(osd + 129);
    const auto *osd_131 = buffer.data(osd + 131);
    const auto *osd_132 = buffer.data(osd + 132);
    const auto *osd_135 = buffer.data(osd + 135);
    const auto *osd_136 = buffer.data(osd + 136);
    const auto *osd_137 = buffer.data(osd + 137);
    const auto *osd_138 = buffer.data(osd + 138);
    const auto *osd_141 = buffer.data(osd + 141);
    const auto *osd_142 = buffer.data(osd + 142);
    const auto *osd_143 = buffer.data(osd + 143);
    const auto *osd_144 = buffer.data(osd + 144);
    const auto *osd_147 = buffer.data(osd + 147);
    const auto *osd_148 = buffer.data(osd + 148);
    const auto *osd_149 = buffer.data(osd + 149);
    const auto *osd_150 = buffer.data(osd + 150);
    const auto *osd_153 = buffer.data(osd + 153);
    const auto *osd_154 = buffer.data(osd + 154);
    const auto *osd_155 = buffer.data(osd + 155);

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pc_x, pc_y, nsd_57, nsd_81, nsd_82, \
                         nsd_83, osp0_40, osp1_40, osd_81, osd_82, \
                         osd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_13 * nsd_81[k]
                   + f_3 * pc_x[k] * osd_81[k];

        t_134[k] = f_13 * nsd_82[k]
                   + f_3 * pc_x[k] * osd_82[k];

        t_135[k] = f_13 * nsd_83[k]
                   + f_3 * pc_x[k] * osd_83[k];

        t_136[k] = f_5 * nsd_57[k]
                   + f_1 * osp0_40[k]
                   - f_2 * osp1_40[k]
                   + f_3 * pc_y[k] * osd_81[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pa_y, pc_y, pc_z, nsf0_99, nsd_51, nsd_59, \
                         nsf1_99, osd_81, osd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_12 * nsd_51[k]
                   + f_3 * pc_z[k] * osd_81[k];

        t_138[k] = f_5 * nsd_59[k]
                   + f_3 * pc_y[k] * osd_83[k];

        t_139[k] = pa_y[k] * nsf0_99[k]
                   - f_4 * pc_y[k] * nsf1_99[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, pc_x, pc_y, pc_z, nsd_54, nsd_84, \
                         nsd_87, osp0_42, osp1_42, osd_84, osd_86, \
                         osd_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_13 * nsd_84[k]
                   + f_1 * osp0_42[k]
                   - f_2 * osp1_42[k]
                   + f_3 * pc_x[k] * osd_84[k];

        t_141[k] = f_3 * pc_y[k] * osd_84[k];

        t_142[k] = f_14 * nsd_54[k]
                   + f_3 * pc_z[k] * osd_84[k];

        t_143[k] = f_13 * nsd_87[k]
                   + f_3 * pc_x[k] * osd_87[k];

        t_144[k] = f_3 * pc_y[k] * osd_86[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, pc_x, pc_y, nsd_89, osp0_43, osp0_44, \
                         osp1_43, osp1_44, osd_87, osd_88, osd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_13 * nsd_89[k]
                   + f_3 * pc_x[k] * osd_89[k];

        t_146[k] = f_1 * osp0_43[k]
                   - f_2 * osp1_43[k]
                   + f_3 * pc_y[k] * osd_87[k];

        t_147[k] = f_7 * osp0_44[k]
                   - f_8 * osp1_44[k]
                   + f_3 * pc_y[k] * osd_88[k];

        t_148[k] = f_3 * pc_y[k] * osd_89[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, pc_x, pc_y, pc_z, nsd_59, nsd_60, nsd_90, \
                         osp0_44, osp0_45, osp1_44, osp1_45, osd_89, \
                         osd_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_14 * nsd_59[k]
                   + f_1 * osp0_44[k]
                   - f_2 * osp1_44[k]
                   + f_3 * pc_z[k] * osd_89[k];

        t_150[k] = f_15 * nsd_90[k]
                   + f_1 * osp0_45[k]
                   - f_2 * osp1_45[k]
                   + f_3 * pc_x[k] * osd_90[k];

        t_151[k] = f_16 * nsd_60[k]
                   + f_3 * pc_y[k] * osd_90[k];

        t_152[k] = f_3 * pc_z[k] * osd_90[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, t_157, pc_x, pc_y, pc_z, nsd_63, nsd_93, \
                         nsd_95, osp0_46, osp1_46, osd_91, osd_93, \
                         osd_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_15 * nsd_93[k]
                   + f_3 * pc_x[k] * osd_93[k];

        t_154[k] = f_3 * pc_z[k] * osd_91[k];

        t_155[k] = f_15 * nsd_95[k]
                   + f_3 * pc_x[k] * osd_95[k];

        t_156[k] = f_16 * nsd_63[k]
                   + f_1 * osp0_46[k]
                   - f_2 * osp1_46[k]
                   + f_3 * pc_y[k] * osd_93[k];

        t_157[k] = f_3 * pc_z[k] * osd_93[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, t_161, pa_z, pc_y, pc_z, nsf0_100, nsd_65, \
                         nsd_66, nsf1_100, osp0_47, osp1_47, osd_95, \
                         osd_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_16 * nsd_65[k]
                   + f_3 * pc_y[k] * osd_95[k];

        t_159[k] = f_1 * osp0_47[k]
                   - f_2 * osp1_47[k]
                   + f_3 * pc_z[k] * osd_95[k];

        t_160[k] = pa_z[k] * nsf0_100[k]
                   - f_4 * pc_z[k] * nsf1_100[k];

        t_161[k] = f_14 * nsd_66[k]
                   + f_3 * pc_y[k] * osd_96[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pc_x, pc_z, nsd_60, nsd_99, nsd_100, \
                         nsd_101, osd_96, osd_99, osd_100, osd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_5 * nsd_60[k]
                   + f_3 * pc_z[k] * osd_96[k];

        t_163[k] = f_15 * nsd_99[k]
                   + f_3 * pc_x[k] * osd_99[k];

        t_164[k] = f_15 * nsd_100[k]
                   + f_3 * pc_x[k] * osd_100[k];

        t_165[k] = f_15 * nsd_101[k]
                   + f_3 * pc_x[k] * osd_101[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pa_z, pc_y, pc_z, nsf0_106, nsd_63, \
                         nsd_65, nsd_71, nsf1_106, osp0_50, osp1_50, osd_99, \
                         osd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = pa_z[k] * nsf0_106[k]
                   - f_4 * pc_z[k] * nsf1_106[k];

        t_167[k] = f_5 * nsd_63[k]
                   + f_3 * pc_z[k] * osd_99[k];

        t_168[k] = f_14 * nsd_71[k]
                   + f_3 * pc_y[k] * osd_101[k];

        t_169[k] = f_5 * nsd_65[k]
                   + f_1 * osp0_50[k]
                   - f_2 * osp1_50[k]
                   + f_3 * pc_z[k] * osd_101[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pc_x, pc_y, pc_z, nsd_66, nsd_72, \
                         nsd_102, nsd_105, osp0_51, osp1_51, osd_102, \
                         osd_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_15 * nsd_102[k]
                   + f_1 * osp0_51[k]
                   - f_2 * osp1_51[k]
                   + f_3 * pc_x[k] * osd_102[k];

        t_171[k] = f_12 * nsd_72[k]
                   + f_3 * pc_y[k] * osd_102[k];

        t_172[k] = f_10 * nsd_66[k]
                   + f_3 * pc_z[k] * osd_102[k];

        t_173[k] = f_15 * nsd_105[k]
                   + f_3 * pc_x[k] * osd_105[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, t_177, pc_x, pc_y, pc_z, nsd_69, nsd_75, \
                         nsd_106, nsd_107, osp0_52, osp1_52, osd_105, osd_106, \
                         osd_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_15 * nsd_106[k]
                   + f_3 * pc_x[k] * osd_106[k];

        t_175[k] = f_15 * nsd_107[k]
                   + f_3 * pc_x[k] * osd_107[k];

        t_176[k] = f_12 * nsd_75[k]
                   + f_1 * osp0_52[k]
                   - f_2 * osp1_52[k]
                   + f_3 * pc_y[k] * osd_105[k];

        t_177[k] = f_10 * nsd_69[k]
                   + f_3 * pc_z[k] * osd_105[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pc_x, pc_y, pc_z, nsd_71, nsd_77, nsd_108, \
                         osp0_53, osp0_54, osp1_53, osp1_54, osd_107, \
                         osd_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_12 * nsd_77[k]
                   + f_3 * pc_y[k] * osd_107[k];

        t_179[k] = f_10 * nsd_71[k]
                   + f_1 * osp0_53[k]
                   - f_2 * osp1_53[k]
                   + f_3 * pc_z[k] * osd_107[k];

        t_180[k] = f_15 * nsd_108[k]
                   + f_1 * osp0_54[k]
                   - f_2 * osp1_54[k]
                   + f_3 * pc_x[k] * osd_108[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pc_x, pc_y, pc_z, nsd_72, nsd_78, \
                         nsd_111, nsd_112, osd_108, osd_111, osd_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_10 * nsd_78[k]
                   + f_3 * pc_y[k] * osd_108[k];

        t_182[k] = f_12 * nsd_72[k]
                   + f_3 * pc_z[k] * osd_108[k];

        t_183[k] = f_15 * nsd_111[k]
                   + f_3 * pc_x[k] * osd_111[k];

        t_184[k] = f_15 * nsd_112[k]
                   + f_3 * pc_x[k] * osd_112[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, pc_y, pc_z, nsd_75, nsd_81, nsd_83, \
                         nsd_113, osp0_55, osp1_55, osd_111, osd_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_15 * nsd_113[k]
                   + f_3 * pc_x[k] * osd_113[k];

        t_186[k] = f_10 * nsd_81[k]
                   + f_1 * osp0_55[k]
                   - f_2 * osp1_55[k]
                   + f_3 * pc_y[k] * osd_111[k];

        t_187[k] = f_12 * nsd_75[k]
                   + f_3 * pc_z[k] * osd_111[k];

        t_188[k] = f_10 * nsd_83[k]
                   + f_3 * pc_y[k] * osd_113[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, t_192, pa_y, pc_y, pc_z, nsf0_140, nsd_77, \
                         nsd_78, nsd_84, nsf1_140, osp0_56, osp1_56, osd_113, \
                         osd_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_12 * nsd_77[k]
                   + f_1 * osp0_56[k]
                   - f_2 * osp1_56[k]
                   + f_3 * pc_z[k] * osd_113[k];

        t_190[k] = pa_y[k] * nsf0_140[k]
                   - f_4 * pc_y[k] * nsf1_140[k];

        t_191[k] = f_5 * nsd_84[k]
                   + f_3 * pc_y[k] * osd_114[k];

        t_192[k] = f_14 * nsd_78[k]
                   + f_3 * pc_z[k] * osd_114[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pc_x, pc_y, nsd_87, nsd_117, nsd_118, \
                         nsd_119, osp0_58, osp1_58, osd_117, osd_118, \
                         osd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_15 * nsd_117[k]
                   + f_3 * pc_x[k] * osd_117[k];

        t_194[k] = f_15 * nsd_118[k]
                   + f_3 * pc_x[k] * osd_118[k];

        t_195[k] = f_15 * nsd_119[k]
                   + f_3 * pc_x[k] * osd_119[k];

        t_196[k] = f_5 * nsd_87[k]
                   + f_1 * osp0_58[k]
                   - f_2 * osp1_58[k]
                   + f_3 * pc_y[k] * osd_117[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, pa_y, pc_y, pc_z, nsf0_149, nsd_81, nsd_89, \
                         nsf1_149, osd_117, osd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_14 * nsd_81[k]
                   + f_3 * pc_z[k] * osd_117[k];

        t_198[k] = f_5 * nsd_89[k]
                   + f_3 * pc_y[k] * osd_119[k];

        t_199[k] = pa_y[k] * nsf0_149[k]
                   - f_4 * pc_y[k] * nsf1_149[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, t_204, pc_x, pc_y, pc_z, nsd_84, nsd_120, \
                         nsd_123, osp0_60, osp1_60, osd_120, osd_122, \
                         osd_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_15 * nsd_120[k]
                   + f_1 * osp0_60[k]
                   - f_2 * osp1_60[k]
                   + f_3 * pc_x[k] * osd_120[k];

        t_201[k] = f_3 * pc_y[k] * osd_120[k];

        t_202[k] = f_16 * nsd_84[k]
                   + f_3 * pc_z[k] * osd_120[k];

        t_203[k] = f_15 * nsd_123[k]
                   + f_3 * pc_x[k] * osd_123[k];

        t_204[k] = f_3 * pc_y[k] * osd_122[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, pc_x, pc_y, nsd_125, osp0_61, osp0_62, \
                         osp1_61, osp1_62, osd_123, osd_124, osd_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_15 * nsd_125[k]
                   + f_3 * pc_x[k] * osd_125[k];

        t_206[k] = f_1 * osp0_61[k]
                   - f_2 * osp1_61[k]
                   + f_3 * pc_y[k] * osd_123[k];

        t_207[k] = f_7 * osp0_62[k]
                   - f_8 * osp1_62[k]
                   + f_3 * pc_y[k] * osd_124[k];

        t_208[k] = f_3 * pc_y[k] * osd_125[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, pc_x, pc_y, pc_z, nsd_89, nsd_90, \
                         nsd_126, osp0_62, osp0_63, osp1_62, osp1_63, osd_125, \
                         osd_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = f_16 * nsd_89[k]
                   + f_1 * osp0_62[k]
                   - f_2 * osp1_62[k]
                   + f_3 * pc_z[k] * osd_125[k];

        t_210[k] = f_16 * nsd_126[k]
                   + f_1 * osp0_63[k]
                   - f_2 * osp1_63[k]
                   + f_3 * pc_x[k] * osd_126[k];

        t_211[k] = f_15 * nsd_90[k]
                   + f_3 * pc_y[k] * osd_126[k];

        t_212[k] = f_3 * pc_z[k] * osd_126[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, t_217, pc_x, pc_y, pc_z, nsd_93, nsd_129, \
                         nsd_131, osp0_64, osp1_64, osd_127, osd_129, \
                         osd_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_16 * nsd_129[k]
                   + f_3 * pc_x[k] * osd_129[k];

        t_214[k] = f_3 * pc_z[k] * osd_127[k];

        t_215[k] = f_16 * nsd_131[k]
                   + f_3 * pc_x[k] * osd_131[k];

        t_216[k] = f_15 * nsd_93[k]
                   + f_1 * osp0_64[k]
                   - f_2 * osp1_64[k]
                   + f_3 * pc_y[k] * osd_129[k];

        t_217[k] = f_3 * pc_z[k] * osd_129[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, pa_z, pc_y, pc_z, nsf0_150, nsd_95, \
                         nsd_96, nsf1_150, osp0_65, osp1_65, osd_131, \
                         osd_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_15 * nsd_95[k]
                   + f_3 * pc_y[k] * osd_131[k];

        t_219[k] = f_1 * osp0_65[k]
                   - f_2 * osp1_65[k]
                   + f_3 * pc_z[k] * osd_131[k];

        t_220[k] = pa_z[k] * nsf0_150[k]
                   - f_4 * pc_z[k] * nsf1_150[k];

        t_221[k] = f_16 * nsd_96[k]
                   + f_3 * pc_y[k] * osd_132[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pc_x, pc_z, nsd_90, nsd_135, nsd_136, \
                         nsd_137, osd_132, osd_135, osd_136, osd_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_5 * nsd_90[k]
                   + f_3 * pc_z[k] * osd_132[k];

        t_223[k] = f_16 * nsd_135[k]
                   + f_3 * pc_x[k] * osd_135[k];

        t_224[k] = f_16 * nsd_136[k]
                   + f_3 * pc_x[k] * osd_136[k];

        t_225[k] = f_16 * nsd_137[k]
                   + f_3 * pc_x[k] * osd_137[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, pa_z, pc_y, pc_z, nsf0_156, nsd_93, \
                         nsd_95, nsd_101, nsf1_156, osp0_68, osp1_68, osd_135, \
                         osd_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = pa_z[k] * nsf0_156[k]
                   - f_4 * pc_z[k] * nsf1_156[k];

        t_227[k] = f_5 * nsd_93[k]
                   + f_3 * pc_z[k] * osd_135[k];

        t_228[k] = f_16 * nsd_101[k]
                   + f_3 * pc_y[k] * osd_137[k];

        t_229[k] = f_5 * nsd_95[k]
                   + f_1 * osp0_68[k]
                   - f_2 * osp1_68[k]
                   + f_3 * pc_z[k] * osd_137[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pc_x, pc_y, pc_z, nsd_96, nsd_102, \
                         nsd_138, nsd_141, osp0_69, osp1_69, osd_138, \
                         osd_141 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_16 * nsd_138[k]
                   + f_1 * osp0_69[k]
                   - f_2 * osp1_69[k]
                   + f_3 * pc_x[k] * osd_138[k];

        t_231[k] = f_14 * nsd_102[k]
                   + f_3 * pc_y[k] * osd_138[k];

        t_232[k] = f_10 * nsd_96[k]
                   + f_3 * pc_z[k] * osd_138[k];

        t_233[k] = f_16 * nsd_141[k]
                   + f_3 * pc_x[k] * osd_141[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, pc_x, pc_y, pc_z, nsd_99, nsd_105, \
                         nsd_142, nsd_143, osp0_70, osp1_70, osd_141, osd_142, \
                         osd_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_16 * nsd_142[k]
                   + f_3 * pc_x[k] * osd_142[k];

        t_235[k] = f_16 * nsd_143[k]
                   + f_3 * pc_x[k] * osd_143[k];

        t_236[k] = f_14 * nsd_105[k]
                   + f_1 * osp0_70[k]
                   - f_2 * osp1_70[k]
                   + f_3 * pc_y[k] * osd_141[k];

        t_237[k] = f_10 * nsd_99[k]
                   + f_3 * pc_z[k] * osd_141[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, pc_x, pc_y, pc_z, nsd_101, nsd_107, nsd_144, \
                         osp0_71, osp0_72, osp1_71, osp1_72, osd_143, \
                         osd_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_14 * nsd_107[k]
                   + f_3 * pc_y[k] * osd_143[k];

        t_239[k] = f_10 * nsd_101[k]
                   + f_1 * osp0_71[k]
                   - f_2 * osp1_71[k]
                   + f_3 * pc_z[k] * osd_143[k];

        t_240[k] = f_16 * nsd_144[k]
                   + f_1 * osp0_72[k]
                   - f_2 * osp1_72[k]
                   + f_3 * pc_x[k] * osd_144[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, pc_x, pc_y, pc_z, nsd_102, nsd_108, \
                         nsd_147, nsd_148, osd_144, osd_147, osd_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = f_12 * nsd_108[k]
                   + f_3 * pc_y[k] * osd_144[k];

        t_242[k] = f_12 * nsd_102[k]
                   + f_3 * pc_z[k] * osd_144[k];

        t_243[k] = f_16 * nsd_147[k]
                   + f_3 * pc_x[k] * osd_147[k];

        t_244[k] = f_16 * nsd_148[k]
                   + f_3 * pc_x[k] * osd_148[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, pc_x, pc_y, pc_z, nsd_105, nsd_111, \
                         nsd_113, nsd_149, osp0_73, osp1_73, osd_147, \
                         osd_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_16 * nsd_149[k]
                   + f_3 * pc_x[k] * osd_149[k];

        t_246[k] = f_12 * nsd_111[k]
                   + f_1 * osp0_73[k]
                   - f_2 * osp1_73[k]
                   + f_3 * pc_y[k] * osd_147[k];

        t_247[k] = f_12 * nsd_105[k]
                   + f_3 * pc_z[k] * osd_147[k];

        t_248[k] = f_12 * nsd_113[k]
                   + f_3 * pc_y[k] * osd_149[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, pc_x, pc_y, pc_z, nsd_107, nsd_114, nsd_150, \
                         osp0_74, osp0_75, osp1_74, osp1_75, osd_149, \
                         osd_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_12 * nsd_107[k]
                   + f_1 * osp0_74[k]
                   - f_2 * osp1_74[k]
                   + f_3 * pc_z[k] * osd_149[k];

        t_250[k] = f_16 * nsd_150[k]
                   + f_1 * osp0_75[k]
                   - f_2 * osp1_75[k]
                   + f_3 * pc_x[k] * osd_150[k];

        t_251[k] = f_10 * nsd_114[k]
                   + f_3 * pc_y[k] * osd_150[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pc_x, pc_z, nsd_108, nsd_153, nsd_154, \
                         nsd_155, osd_150, osd_153, osd_154, osd_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_14 * nsd_108[k]
                   + f_3 * pc_z[k] * osd_150[k];

        t_253[k] = f_16 * nsd_153[k]
                   + f_3 * pc_x[k] * osd_153[k];

        t_254[k] = f_16 * nsd_154[k]
                   + f_3 * pc_x[k] * osd_154[k];

        t_255[k] = f_16 * nsd_155[k]
                   + f_3 * pc_x[k] * osd_155[k];
    }
}

static auto
compute_prim_osf_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsf0,
                                                          const size_t nsd, const size_t nsf1,
                                                          const size_t osp0, const size_t osp1,
                                                          const size_t osd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_10 = 1.0 / q;
    const auto f_11 = 4.0 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 3.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.5 / q;

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
    auto *t_373 = buffer.data(target + 373);
    auto *t_374 = buffer.data(target + 374);
    auto *t_375 = buffer.data(target + 375);
    auto *t_376 = buffer.data(target + 376);
    auto *t_377 = buffer.data(target + 377);
    auto *t_378 = buffer.data(target + 378);
    auto *t_379 = buffer.data(target + 379);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsf0_200 = buffer.data(nsf0 + 200);
    const auto *nsf0_209 = buffer.data(nsf0 + 209);
    const auto *nsf0_210 = buffer.data(nsf0 + 210);
    const auto *nsf0_216 = buffer.data(nsf0 + 216);
    const auto *nsf0_270 = buffer.data(nsf0 + 270);
    const auto *nsf0_279 = buffer.data(nsf0 + 279);
    const auto *nsf0_280 = buffer.data(nsf0 + 280);
    const auto *nsf0_286 = buffer.data(nsf0 + 286);

    const auto *nsd_111 = buffer.data(nsd + 111);
    const auto *nsd_113 = buffer.data(nsd + 113);
    const auto *nsd_114 = buffer.data(nsd + 114);
    const auto *nsd_117 = buffer.data(nsd + 117);
    const auto *nsd_119 = buffer.data(nsd + 119);
    const auto *nsd_120 = buffer.data(nsd + 120);
    const auto *nsd_123 = buffer.data(nsd + 123);
    const auto *nsd_125 = buffer.data(nsd + 125);
    const auto *nsd_126 = buffer.data(nsd + 126);
    const auto *nsd_129 = buffer.data(nsd + 129);
    const auto *nsd_131 = buffer.data(nsd + 131);
    const auto *nsd_132 = buffer.data(nsd + 132);
    const auto *nsd_135 = buffer.data(nsd + 135);
    const auto *nsd_137 = buffer.data(nsd + 137);
    const auto *nsd_138 = buffer.data(nsd + 138);
    const auto *nsd_141 = buffer.data(nsd + 141);
    const auto *nsd_143 = buffer.data(nsd + 143);
    const auto *nsd_144 = buffer.data(nsd + 144);
    const auto *nsd_147 = buffer.data(nsd + 147);
    const auto *nsd_149 = buffer.data(nsd + 149);
    const auto *nsd_150 = buffer.data(nsd + 150);
    const auto *nsd_153 = buffer.data(nsd + 153);
    const auto *nsd_155 = buffer.data(nsd + 155);
    const auto *nsd_156 = buffer.data(nsd + 156);
    const auto *nsd_159 = buffer.data(nsd + 159);
    const auto *nsd_160 = buffer.data(nsd + 160);
    const auto *nsd_161 = buffer.data(nsd + 161);
    const auto *nsd_162 = buffer.data(nsd + 162);
    const auto *nsd_165 = buffer.data(nsd + 165);
    const auto *nsd_167 = buffer.data(nsd + 167);
    const auto *nsd_168 = buffer.data(nsd + 168);
    const auto *nsd_171 = buffer.data(nsd + 171);
    const auto *nsd_173 = buffer.data(nsd + 173);
    const auto *nsd_174 = buffer.data(nsd + 174);
    const auto *nsd_177 = buffer.data(nsd + 177);
    const auto *nsd_178 = buffer.data(nsd + 178);
    const auto *nsd_179 = buffer.data(nsd + 179);
    const auto *nsd_180 = buffer.data(nsd + 180);
    const auto *nsd_183 = buffer.data(nsd + 183);
    const auto *nsd_184 = buffer.data(nsd + 184);
    const auto *nsd_185 = buffer.data(nsd + 185);
    const auto *nsd_186 = buffer.data(nsd + 186);
    const auto *nsd_189 = buffer.data(nsd + 189);
    const auto *nsd_190 = buffer.data(nsd + 190);
    const auto *nsd_191 = buffer.data(nsd + 191);
    const auto *nsd_192 = buffer.data(nsd + 192);
    const auto *nsd_195 = buffer.data(nsd + 195);
    const auto *nsd_196 = buffer.data(nsd + 196);
    const auto *nsd_197 = buffer.data(nsd + 197);
    const auto *nsd_198 = buffer.data(nsd + 198);
    const auto *nsd_201 = buffer.data(nsd + 201);
    const auto *nsd_202 = buffer.data(nsd + 202);
    const auto *nsd_203 = buffer.data(nsd + 203);
    const auto *nsd_207 = buffer.data(nsd + 207);
    const auto *nsd_208 = buffer.data(nsd + 208);
    const auto *nsd_209 = buffer.data(nsd + 209);
    const auto *nsd_210 = buffer.data(nsd + 210);
    const auto *nsd_213 = buffer.data(nsd + 213);
    const auto *nsd_215 = buffer.data(nsd + 215);
    const auto *nsd_216 = buffer.data(nsd + 216);
    const auto *nsd_219 = buffer.data(nsd + 219);
    const auto *nsd_221 = buffer.data(nsd + 221);
    const auto *nsd_225 = buffer.data(nsd + 225);
    const auto *nsd_226 = buffer.data(nsd + 226);
    const auto *nsd_227 = buffer.data(nsd + 227);

    const auto *nsf1_200 = buffer.data(nsf1 + 200);
    const auto *nsf1_209 = buffer.data(nsf1 + 209);
    const auto *nsf1_210 = buffer.data(nsf1 + 210);
    const auto *nsf1_216 = buffer.data(nsf1 + 216);
    const auto *nsf1_270 = buffer.data(nsf1 + 270);
    const auto *nsf1_279 = buffer.data(nsf1 + 279);
    const auto *nsf1_280 = buffer.data(nsf1 + 280);
    const auto *nsf1_286 = buffer.data(nsf1 + 286);

    const auto *osp0_76 = buffer.data(osp0 + 76);
    const auto *osp0_77 = buffer.data(osp0 + 77);
    const auto *osp0_79 = buffer.data(osp0 + 79);
    const auto *osp0_81 = buffer.data(osp0 + 81);
    const auto *osp0_82 = buffer.data(osp0 + 82);
    const auto *osp0_83 = buffer.data(osp0 + 83);
    const auto *osp0_84 = buffer.data(osp0 + 84);
    const auto *osp0_85 = buffer.data(osp0 + 85);
    const auto *osp0_86 = buffer.data(osp0 + 86);
    const auto *osp0_89 = buffer.data(osp0 + 89);
    const auto *osp0_90 = buffer.data(osp0 + 90);
    const auto *osp0_91 = buffer.data(osp0 + 91);
    const auto *osp0_92 = buffer.data(osp0 + 92);
    const auto *osp0_93 = buffer.data(osp0 + 93);
    const auto *osp0_94 = buffer.data(osp0 + 94);
    const auto *osp0_95 = buffer.data(osp0 + 95);
    const auto *osp0_96 = buffer.data(osp0 + 96);
    const auto *osp0_97 = buffer.data(osp0 + 97);
    const auto *osp0_98 = buffer.data(osp0 + 98);
    const auto *osp0_99 = buffer.data(osp0 + 99);
    const auto *osp0_100 = buffer.data(osp0 + 100);
    const auto *osp0_101 = buffer.data(osp0 + 101);
    const auto *osp0_103 = buffer.data(osp0 + 103);
    const auto *osp0_105 = buffer.data(osp0 + 105);
    const auto *osp0_106 = buffer.data(osp0 + 106);
    const auto *osp0_107 = buffer.data(osp0 + 107);
    const auto *osp0_108 = buffer.data(osp0 + 108);
    const auto *osp0_109 = buffer.data(osp0 + 109);
    const auto *osp0_110 = buffer.data(osp0 + 110);
    const auto *osp0_113 = buffer.data(osp0 + 113);

    const auto *osp1_76 = buffer.data(osp1 + 76);
    const auto *osp1_77 = buffer.data(osp1 + 77);
    const auto *osp1_79 = buffer.data(osp1 + 79);
    const auto *osp1_81 = buffer.data(osp1 + 81);
    const auto *osp1_82 = buffer.data(osp1 + 82);
    const auto *osp1_83 = buffer.data(osp1 + 83);
    const auto *osp1_84 = buffer.data(osp1 + 84);
    const auto *osp1_85 = buffer.data(osp1 + 85);
    const auto *osp1_86 = buffer.data(osp1 + 86);
    const auto *osp1_89 = buffer.data(osp1 + 89);
    const auto *osp1_90 = buffer.data(osp1 + 90);
    const auto *osp1_91 = buffer.data(osp1 + 91);
    const auto *osp1_92 = buffer.data(osp1 + 92);
    const auto *osp1_93 = buffer.data(osp1 + 93);
    const auto *osp1_94 = buffer.data(osp1 + 94);
    const auto *osp1_95 = buffer.data(osp1 + 95);
    const auto *osp1_96 = buffer.data(osp1 + 96);
    const auto *osp1_97 = buffer.data(osp1 + 97);
    const auto *osp1_98 = buffer.data(osp1 + 98);
    const auto *osp1_99 = buffer.data(osp1 + 99);
    const auto *osp1_100 = buffer.data(osp1 + 100);
    const auto *osp1_101 = buffer.data(osp1 + 101);
    const auto *osp1_103 = buffer.data(osp1 + 103);
    const auto *osp1_105 = buffer.data(osp1 + 105);
    const auto *osp1_106 = buffer.data(osp1 + 106);
    const auto *osp1_107 = buffer.data(osp1 + 107);
    const auto *osp1_108 = buffer.data(osp1 + 108);
    const auto *osp1_109 = buffer.data(osp1 + 109);
    const auto *osp1_110 = buffer.data(osp1 + 110);
    const auto *osp1_113 = buffer.data(osp1 + 113);

    const auto *osd_153 = buffer.data(osd + 153);
    const auto *osd_155 = buffer.data(osd + 155);
    const auto *osd_156 = buffer.data(osd + 156);
    const auto *osd_159 = buffer.data(osd + 159);
    const auto *osd_160 = buffer.data(osd + 160);
    const auto *osd_161 = buffer.data(osd + 161);
    const auto *osd_162 = buffer.data(osd + 162);
    const auto *osd_164 = buffer.data(osd + 164);
    const auto *osd_165 = buffer.data(osd + 165);
    const auto *osd_166 = buffer.data(osd + 166);
    const auto *osd_167 = buffer.data(osd + 167);
    const auto *osd_168 = buffer.data(osd + 168);
    const auto *osd_169 = buffer.data(osd + 169);
    const auto *osd_171 = buffer.data(osd + 171);
    const auto *osd_173 = buffer.data(osd + 173);
    const auto *osd_174 = buffer.data(osd + 174);
    const auto *osd_177 = buffer.data(osd + 177);
    const auto *osd_178 = buffer.data(osd + 178);
    const auto *osd_179 = buffer.data(osd + 179);
    const auto *osd_180 = buffer.data(osd + 180);
    const auto *osd_183 = buffer.data(osd + 183);
    const auto *osd_184 = buffer.data(osd + 184);
    const auto *osd_185 = buffer.data(osd + 185);
    const auto *osd_186 = buffer.data(osd + 186);
    const auto *osd_189 = buffer.data(osd + 189);
    const auto *osd_190 = buffer.data(osd + 190);
    const auto *osd_191 = buffer.data(osd + 191);
    const auto *osd_192 = buffer.data(osd + 192);
    const auto *osd_195 = buffer.data(osd + 195);
    const auto *osd_196 = buffer.data(osd + 196);
    const auto *osd_197 = buffer.data(osd + 197);
    const auto *osd_198 = buffer.data(osd + 198);
    const auto *osd_201 = buffer.data(osd + 201);
    const auto *osd_202 = buffer.data(osd + 202);
    const auto *osd_203 = buffer.data(osd + 203);
    const auto *osd_204 = buffer.data(osd + 204);
    const auto *osd_207 = buffer.data(osd + 207);
    const auto *osd_208 = buffer.data(osd + 208);
    const auto *osd_209 = buffer.data(osd + 209);
    const auto *osd_210 = buffer.data(osd + 210);
    const auto *osd_212 = buffer.data(osd + 212);
    const auto *osd_213 = buffer.data(osd + 213);
    const auto *osd_214 = buffer.data(osd + 214);
    const auto *osd_215 = buffer.data(osd + 215);
    const auto *osd_216 = buffer.data(osd + 216);
    const auto *osd_217 = buffer.data(osd + 217);
    const auto *osd_219 = buffer.data(osd + 219);
    const auto *osd_221 = buffer.data(osd + 221);
    const auto *osd_222 = buffer.data(osd + 222);
    const auto *osd_225 = buffer.data(osd + 225);
    const auto *osd_226 = buffer.data(osd + 226);
    const auto *osd_227 = buffer.data(osd + 227);

#pragma omp simd aligned(t_256, t_257, t_258, t_259, pc_y, pc_z, nsd_111, nsd_113, nsd_117, \
                         nsd_119, osp0_76, osp0_77, osp1_76, osp1_77, osd_153, \
                         osd_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_10 * nsd_117[k]
                   + f_1 * osp0_76[k]
                   - f_2 * osp1_76[k]
                   + f_3 * pc_y[k] * osd_153[k];

        t_257[k] = f_14 * nsd_111[k]
                   + f_3 * pc_z[k] * osd_153[k];

        t_258[k] = f_10 * nsd_119[k]
                   + f_3 * pc_y[k] * osd_155[k];

        t_259[k] = f_14 * nsd_113[k]
                   + f_1 * osp0_77[k]
                   - f_2 * osp1_77[k]
                   + f_3 * pc_z[k] * osd_155[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pa_y, pc_x, pc_y, pc_z, nsf0_200, \
                         nsd_114, nsd_120, nsd_159, nsf1_200, osd_156, \
                         osd_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = pa_y[k] * nsf0_200[k]
                   - f_4 * pc_y[k] * nsf1_200[k];

        t_261[k] = f_5 * nsd_120[k]
                   + f_3 * pc_y[k] * osd_156[k];

        t_262[k] = f_16 * nsd_114[k]
                   + f_3 * pc_z[k] * osd_156[k];

        t_263[k] = f_16 * nsd_159[k]
                   + f_3 * pc_x[k] * osd_159[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pc_x, pc_y, pc_z, nsd_117, nsd_123, \
                         nsd_160, nsd_161, osp0_79, osp1_79, osd_159, osd_160, \
                         osd_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_16 * nsd_160[k]
                   + f_3 * pc_x[k] * osd_160[k];

        t_265[k] = f_16 * nsd_161[k]
                   + f_3 * pc_x[k] * osd_161[k];

        t_266[k] = f_5 * nsd_123[k]
                   + f_1 * osp0_79[k]
                   - f_2 * osp1_79[k]
                   + f_3 * pc_y[k] * osd_159[k];

        t_267[k] = f_16 * nsd_117[k]
                   + f_3 * pc_z[k] * osd_159[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, t_271, pa_y, pc_x, pc_y, nsf0_209, nsd_125, \
                         nsd_162, nsf1_209, osp0_81, osp1_81, osd_161, \
                         osd_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_5 * nsd_125[k]
                   + f_3 * pc_y[k] * osd_161[k];

        t_269[k] = pa_y[k] * nsf0_209[k]
                   - f_4 * pc_y[k] * nsf1_209[k];

        t_270[k] = f_16 * nsd_162[k]
                   + f_1 * osp0_81[k]
                   - f_2 * osp1_81[k]
                   + f_3 * pc_x[k] * osd_162[k];

        t_271[k] = f_3 * pc_y[k] * osd_162[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pc_x, pc_y, pc_z, nsd_120, nsd_165, \
                         nsd_167, osd_162, osd_164, osd_165, osd_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = f_15 * nsd_120[k]
                   + f_3 * pc_z[k] * osd_162[k];

        t_273[k] = f_16 * nsd_165[k]
                   + f_3 * pc_x[k] * osd_165[k];

        t_274[k] = f_3 * pc_y[k] * osd_164[k];

        t_275[k] = f_16 * nsd_167[k]
                   + f_3 * pc_x[k] * osd_167[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_y, pc_z, nsd_125, osp0_82, osp0_83, \
                         osp1_82, osp1_83, osd_165, osd_166, osd_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_1 * osp0_82[k]
                   - f_2 * osp1_82[k]
                   + f_3 * pc_y[k] * osd_165[k];

        t_277[k] = f_7 * osp0_83[k]
                   - f_8 * osp1_83[k]
                   + f_3 * pc_y[k] * osd_166[k];

        t_278[k] = f_3 * pc_y[k] * osd_167[k];

        t_279[k] = f_15 * nsd_125[k]
                   + f_1 * osp0_83[k]
                   - f_2 * osp1_83[k]
                   + f_3 * pc_z[k] * osd_167[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, t_284, pc_x, pc_y, pc_z, nsd_126, \
                         nsd_168, nsd_171, osp0_84, osp1_84, osd_168, osd_169, \
                         osd_171 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_14 * nsd_168[k]
                   + f_1 * osp0_84[k]
                   - f_2 * osp1_84[k]
                   + f_3 * pc_x[k] * osd_168[k];

        t_281[k] = f_13 * nsd_126[k]
                   + f_3 * pc_y[k] * osd_168[k];

        t_282[k] = f_3 * pc_z[k] * osd_168[k];

        t_283[k] = f_14 * nsd_171[k]
                   + f_3 * pc_x[k] * osd_171[k];

        t_284[k] = f_3 * pc_z[k] * osd_169[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, pc_x, pc_y, pc_z, nsd_129, nsd_131, \
                         nsd_173, osp0_85, osp1_85, osd_171, osd_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = f_14 * nsd_173[k]
                   + f_3 * pc_x[k] * osd_173[k];

        t_286[k] = f_13 * nsd_129[k]
                   + f_1 * osp0_85[k]
                   - f_2 * osp1_85[k]
                   + f_3 * pc_y[k] * osd_171[k];

        t_287[k] = f_3 * pc_z[k] * osd_171[k];

        t_288[k] = f_13 * nsd_131[k]
                   + f_3 * pc_y[k] * osd_173[k];
    }

#pragma omp simd aligned(t_289, t_290, t_291, t_292, pa_z, pc_y, pc_z, nsf0_210, nsd_126, \
                         nsd_132, nsf1_210, osp0_86, osp1_86, osd_173, \
                         osd_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_289[k] = f_1 * osp0_86[k]
                   - f_2 * osp1_86[k]
                   + f_3 * pc_z[k] * osd_173[k];

        t_290[k] = pa_z[k] * nsf0_210[k]
                   - f_4 * pc_z[k] * nsf1_210[k];

        t_291[k] = f_15 * nsd_132[k]
                   + f_3 * pc_y[k] * osd_174[k];

        t_292[k] = f_5 * nsd_126[k]
                   + f_3 * pc_z[k] * osd_174[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, pa_z, pc_x, pc_z, nsf0_216, nsd_177, \
                         nsd_178, nsd_179, nsf1_216, osd_177, osd_178, \
                         osd_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = f_14 * nsd_177[k]
                   + f_3 * pc_x[k] * osd_177[k];

        t_294[k] = f_14 * nsd_178[k]
                   + f_3 * pc_x[k] * osd_178[k];

        t_295[k] = f_14 * nsd_179[k]
                   + f_3 * pc_x[k] * osd_179[k];

        t_296[k] = pa_z[k] * nsf0_216[k]
                   - f_4 * pc_z[k] * nsf1_216[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, pc_y, pc_z, nsd_129, nsd_131, nsd_137, osp0_89, \
                         osp1_89, osd_177, osd_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_5 * nsd_129[k]
                   + f_3 * pc_z[k] * osd_177[k];

        t_298[k] = f_15 * nsd_137[k]
                   + f_3 * pc_y[k] * osd_179[k];

        t_299[k] = f_5 * nsd_131[k]
                   + f_1 * osp0_89[k]
                   - f_2 * osp1_89[k]
                   + f_3 * pc_z[k] * osd_179[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, pc_x, pc_y, pc_z, nsd_132, nsd_138, \
                         nsd_180, nsd_183, osp0_90, osp1_90, osd_180, \
                         osd_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_14 * nsd_180[k]
                   + f_1 * osp0_90[k]
                   - f_2 * osp1_90[k]
                   + f_3 * pc_x[k] * osd_180[k];

        t_301[k] = f_16 * nsd_138[k]
                   + f_3 * pc_y[k] * osd_180[k];

        t_302[k] = f_10 * nsd_132[k]
                   + f_3 * pc_z[k] * osd_180[k];

        t_303[k] = f_14 * nsd_183[k]
                   + f_3 * pc_x[k] * osd_183[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, pc_x, pc_y, pc_z, nsd_135, nsd_141, \
                         nsd_184, nsd_185, osp0_91, osp1_91, osd_183, osd_184, \
                         osd_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_14 * nsd_184[k]
                   + f_3 * pc_x[k] * osd_184[k];

        t_305[k] = f_14 * nsd_185[k]
                   + f_3 * pc_x[k] * osd_185[k];

        t_306[k] = f_16 * nsd_141[k]
                   + f_1 * osp0_91[k]
                   - f_2 * osp1_91[k]
                   + f_3 * pc_y[k] * osd_183[k];

        t_307[k] = f_10 * nsd_135[k]
                   + f_3 * pc_z[k] * osd_183[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, pc_x, pc_y, pc_z, nsd_137, nsd_143, nsd_186, \
                         osp0_92, osp0_93, osp1_92, osp1_93, osd_185, \
                         osd_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = f_16 * nsd_143[k]
                   + f_3 * pc_y[k] * osd_185[k];

        t_309[k] = f_10 * nsd_137[k]
                   + f_1 * osp0_92[k]
                   - f_2 * osp1_92[k]
                   + f_3 * pc_z[k] * osd_185[k];

        t_310[k] = f_14 * nsd_186[k]
                   + f_1 * osp0_93[k]
                   - f_2 * osp1_93[k]
                   + f_3 * pc_x[k] * osd_186[k];
    }

#pragma omp simd aligned(t_311, t_312, t_313, t_314, pc_x, pc_y, pc_z, nsd_138, nsd_144, \
                         nsd_189, nsd_190, osd_186, osd_189, osd_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_311[k] = f_14 * nsd_144[k]
                   + f_3 * pc_y[k] * osd_186[k];

        t_312[k] = f_12 * nsd_138[k]
                   + f_3 * pc_z[k] * osd_186[k];

        t_313[k] = f_14 * nsd_189[k]
                   + f_3 * pc_x[k] * osd_189[k];

        t_314[k] = f_14 * nsd_190[k]
                   + f_3 * pc_x[k] * osd_190[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, t_318, pc_x, pc_y, pc_z, nsd_141, nsd_147, \
                         nsd_149, nsd_191, osp0_94, osp1_94, osd_189, \
                         osd_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_14 * nsd_191[k]
                   + f_3 * pc_x[k] * osd_191[k];

        t_316[k] = f_14 * nsd_147[k]
                   + f_1 * osp0_94[k]
                   - f_2 * osp1_94[k]
                   + f_3 * pc_y[k] * osd_189[k];

        t_317[k] = f_12 * nsd_141[k]
                   + f_3 * pc_z[k] * osd_189[k];

        t_318[k] = f_14 * nsd_149[k]
                   + f_3 * pc_y[k] * osd_191[k];
    }

#pragma omp simd aligned(t_319, t_320, t_321, pc_x, pc_y, pc_z, nsd_143, nsd_150, nsd_192, \
                         osp0_95, osp0_96, osp1_95, osp1_96, osd_191, \
                         osd_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_12 * nsd_143[k]
                   + f_1 * osp0_95[k]
                   - f_2 * osp1_95[k]
                   + f_3 * pc_z[k] * osd_191[k];

        t_320[k] = f_14 * nsd_192[k]
                   + f_1 * osp0_96[k]
                   - f_2 * osp1_96[k]
                   + f_3 * pc_x[k] * osd_192[k];

        t_321[k] = f_12 * nsd_150[k]
                   + f_3 * pc_y[k] * osd_192[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, t_325, pc_x, pc_z, nsd_144, nsd_195, nsd_196, \
                         nsd_197, osd_192, osd_195, osd_196, osd_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_14 * nsd_144[k]
                   + f_3 * pc_z[k] * osd_192[k];

        t_323[k] = f_14 * nsd_195[k]
                   + f_3 * pc_x[k] * osd_195[k];

        t_324[k] = f_14 * nsd_196[k]
                   + f_3 * pc_x[k] * osd_196[k];

        t_325[k] = f_14 * nsd_197[k]
                   + f_3 * pc_x[k] * osd_197[k];
    }

#pragma omp simd aligned(t_326, t_327, t_328, t_329, pc_y, pc_z, nsd_147, nsd_149, nsd_153, \
                         nsd_155, osp0_97, osp0_98, osp1_97, osp1_98, osd_195, \
                         osd_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_326[k] = f_12 * nsd_153[k]
                   + f_1 * osp0_97[k]
                   - f_2 * osp1_97[k]
                   + f_3 * pc_y[k] * osd_195[k];

        t_327[k] = f_14 * nsd_147[k]
                   + f_3 * pc_z[k] * osd_195[k];

        t_328[k] = f_12 * nsd_155[k]
                   + f_3 * pc_y[k] * osd_197[k];

        t_329[k] = f_14 * nsd_149[k]
                   + f_1 * osp0_98[k]
                   - f_2 * osp1_98[k]
                   + f_3 * pc_z[k] * osd_197[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, t_333, pc_x, pc_y, pc_z, nsd_150, nsd_156, \
                         nsd_198, nsd_201, osp0_99, osp1_99, osd_198, \
                         osd_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = f_14 * nsd_198[k]
                   + f_1 * osp0_99[k]
                   - f_2 * osp1_99[k]
                   + f_3 * pc_x[k] * osd_198[k];

        t_331[k] = f_10 * nsd_156[k]
                   + f_3 * pc_y[k] * osd_198[k];

        t_332[k] = f_16 * nsd_150[k]
                   + f_3 * pc_z[k] * osd_198[k];

        t_333[k] = f_14 * nsd_201[k]
                   + f_3 * pc_x[k] * osd_201[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, t_337, pc_x, pc_y, pc_z, nsd_153, nsd_159, \
                         nsd_202, nsd_203, osp0_100, osp1_100, osd_201, osd_202, \
                         osd_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_14 * nsd_202[k]
                   + f_3 * pc_x[k] * osd_202[k];

        t_335[k] = f_14 * nsd_203[k]
                   + f_3 * pc_x[k] * osd_203[k];

        t_336[k] = f_10 * nsd_159[k]
                   + f_1 * osp0_100[k]
                   - f_2 * osp1_100[k]
                   + f_3 * pc_y[k] * osd_201[k];

        t_337[k] = f_16 * nsd_153[k]
                   + f_3 * pc_z[k] * osd_201[k];
    }

#pragma omp simd aligned(t_338, t_339, t_340, t_341, pa_y, pc_y, pc_z, nsf0_270, nsd_155, \
                         nsd_161, nsd_162, nsf1_270, osp0_101, osp1_101, osd_203, \
                         osd_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_338[k] = f_10 * nsd_161[k]
                   + f_3 * pc_y[k] * osd_203[k];

        t_339[k] = f_16 * nsd_155[k]
                   + f_1 * osp0_101[k]
                   - f_2 * osp1_101[k]
                   + f_3 * pc_z[k] * osd_203[k];

        t_340[k] = pa_y[k] * nsf0_270[k]
                   - f_4 * pc_y[k] * nsf1_270[k];

        t_341[k] = f_5 * nsd_162[k]
                   + f_3 * pc_y[k] * osd_204[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, t_345, pc_x, pc_z, nsd_156, nsd_207, nsd_208, \
                         nsd_209, osd_204, osd_207, osd_208, osd_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_15 * nsd_156[k]
                   + f_3 * pc_z[k] * osd_204[k];

        t_343[k] = f_14 * nsd_207[k]
                   + f_3 * pc_x[k] * osd_207[k];

        t_344[k] = f_14 * nsd_208[k]
                   + f_3 * pc_x[k] * osd_208[k];

        t_345[k] = f_14 * nsd_209[k]
                   + f_3 * pc_x[k] * osd_209[k];
    }

#pragma omp simd aligned(t_346, t_347, t_348, t_349, pa_y, pc_y, pc_z, nsf0_279, nsd_159, \
                         nsd_165, nsd_167, nsf1_279, osp0_103, osp1_103, osd_207, \
                         osd_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_346[k] = f_5 * nsd_165[k]
                   + f_1 * osp0_103[k]
                   - f_2 * osp1_103[k]
                   + f_3 * pc_y[k] * osd_207[k];

        t_347[k] = f_15 * nsd_159[k]
                   + f_3 * pc_z[k] * osd_207[k];

        t_348[k] = f_5 * nsd_167[k]
                   + f_3 * pc_y[k] * osd_209[k];

        t_349[k] = pa_y[k] * nsf0_279[k]
                   - f_4 * pc_y[k] * nsf1_279[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, t_354, pc_x, pc_y, pc_z, nsd_162, \
                         nsd_210, nsd_213, osp0_105, osp1_105, osd_210, osd_212, \
                         osd_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_14 * nsd_210[k]
                   + f_1 * osp0_105[k]
                   - f_2 * osp1_105[k]
                   + f_3 * pc_x[k] * osd_210[k];

        t_351[k] = f_3 * pc_y[k] * osd_210[k];

        t_352[k] = f_13 * nsd_162[k]
                   + f_3 * pc_z[k] * osd_210[k];

        t_353[k] = f_14 * nsd_213[k]
                   + f_3 * pc_x[k] * osd_213[k];

        t_354[k] = f_3 * pc_y[k] * osd_212[k];
    }

#pragma omp simd aligned(t_355, t_356, t_357, t_358, pc_x, pc_y, nsd_215, osp0_106, osp0_107, \
                         osp1_106, osp1_107, osd_213, osd_214, \
                         osd_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_355[k] = f_14 * nsd_215[k]
                   + f_3 * pc_x[k] * osd_215[k];

        t_356[k] = f_1 * osp0_106[k]
                   - f_2 * osp1_106[k]
                   + f_3 * pc_y[k] * osd_213[k];

        t_357[k] = f_7 * osp0_107[k]
                   - f_8 * osp1_107[k]
                   + f_3 * pc_y[k] * osd_214[k];

        t_358[k] = f_3 * pc_y[k] * osd_215[k];
    }

#pragma omp simd aligned(t_359, t_360, t_361, t_362, pc_x, pc_y, pc_z, nsd_167, nsd_168, \
                         nsd_216, osp0_107, osp0_108, osp1_107, osp1_108, osd_215, \
                         osd_216 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = f_13 * nsd_167[k]
                   + f_1 * osp0_107[k]
                   - f_2 * osp1_107[k]
                   + f_3 * pc_z[k] * osd_215[k];

        t_360[k] = f_12 * nsd_216[k]
                   + f_1 * osp0_108[k]
                   - f_2 * osp1_108[k]
                   + f_3 * pc_x[k] * osd_216[k];

        t_361[k] = f_11 * nsd_168[k]
                   + f_3 * pc_y[k] * osd_216[k];

        t_362[k] = f_3 * pc_z[k] * osd_216[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, t_366, t_367, pc_x, pc_y, pc_z, nsd_171, \
                         nsd_219, nsd_221, osp0_109, osp1_109, osd_217, osd_219, \
                         osd_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = f_12 * nsd_219[k]
                   + f_3 * pc_x[k] * osd_219[k];

        t_364[k] = f_3 * pc_z[k] * osd_217[k];

        t_365[k] = f_12 * nsd_221[k]
                   + f_3 * pc_x[k] * osd_221[k];

        t_366[k] = f_11 * nsd_171[k]
                   + f_1 * osp0_109[k]
                   - f_2 * osp1_109[k]
                   + f_3 * pc_y[k] * osd_219[k];

        t_367[k] = f_3 * pc_z[k] * osd_219[k];
    }

#pragma omp simd aligned(t_368, t_369, t_370, t_371, pa_z, pc_y, pc_z, nsf0_280, nsd_173, \
                         nsd_174, nsf1_280, osp0_110, osp1_110, osd_221, \
                         osd_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_368[k] = f_11 * nsd_173[k]
                   + f_3 * pc_y[k] * osd_221[k];

        t_369[k] = f_1 * osp0_110[k]
                   - f_2 * osp1_110[k]
                   + f_3 * pc_z[k] * osd_221[k];

        t_370[k] = pa_z[k] * nsf0_280[k]
                   - f_4 * pc_z[k] * nsf1_280[k];

        t_371[k] = f_13 * nsd_174[k]
                   + f_3 * pc_y[k] * osd_222[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, t_375, pc_x, pc_z, nsd_168, nsd_225, nsd_226, \
                         nsd_227, osd_222, osd_225, osd_226, osd_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = f_5 * nsd_168[k]
                   + f_3 * pc_z[k] * osd_222[k];

        t_373[k] = f_12 * nsd_225[k]
                   + f_3 * pc_x[k] * osd_225[k];

        t_374[k] = f_12 * nsd_226[k]
                   + f_3 * pc_x[k] * osd_226[k];

        t_375[k] = f_12 * nsd_227[k]
                   + f_3 * pc_x[k] * osd_227[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, t_379, pa_z, pc_y, pc_z, nsf0_286, nsd_171, \
                         nsd_173, nsd_179, nsf1_286, osp0_113, osp1_113, osd_225, \
                         osd_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = pa_z[k] * nsf0_286[k]
                   - f_4 * pc_z[k] * nsf1_286[k];

        t_377[k] = f_5 * nsd_171[k]
                   + f_3 * pc_z[k] * osd_225[k];

        t_378[k] = f_13 * nsd_179[k]
                   + f_3 * pc_y[k] * osd_227[k];

        t_379[k] = f_5 * nsd_173[k]
                   + f_1 * osp0_113[k]
                   - f_2 * osp1_113[k]
                   + f_3 * pc_z[k] * osd_227[k];
    }
}

static auto
compute_prim_osf_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsf0,
                                                          const size_t nsd, const size_t nsf1,
                                                          const size_t osp0, const size_t osp1,
                                                          const size_t osd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_9 = 4.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 4.0 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 3.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.5 / q;

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
    auto *t_492 = buffer.data(target + 492);
    auto *t_493 = buffer.data(target + 493);
    auto *t_494 = buffer.data(target + 494);
    auto *t_495 = buffer.data(target + 495);
    auto *t_496 = buffer.data(target + 496);
    auto *t_497 = buffer.data(target + 497);
    auto *t_498 = buffer.data(target + 498);
    auto *t_499 = buffer.data(target + 499);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsf0_350 = buffer.data(nsf0 + 350);
    const auto *nsf0_359 = buffer.data(nsf0 + 359);
    const auto *nsf0_360 = buffer.data(nsf0 + 360);
    const auto *nsf0_366 = buffer.data(nsf0 + 366);

    const auto *nsd_174 = buffer.data(nsd + 174);
    const auto *nsd_177 = buffer.data(nsd + 177);
    const auto *nsd_179 = buffer.data(nsd + 179);
    const auto *nsd_180 = buffer.data(nsd + 180);
    const auto *nsd_183 = buffer.data(nsd + 183);
    const auto *nsd_185 = buffer.data(nsd + 185);
    const auto *nsd_186 = buffer.data(nsd + 186);
    const auto *nsd_189 = buffer.data(nsd + 189);
    const auto *nsd_191 = buffer.data(nsd + 191);
    const auto *nsd_192 = buffer.data(nsd + 192);
    const auto *nsd_195 = buffer.data(nsd + 195);
    const auto *nsd_197 = buffer.data(nsd + 197);
    const auto *nsd_198 = buffer.data(nsd + 198);
    const auto *nsd_201 = buffer.data(nsd + 201);
    const auto *nsd_203 = buffer.data(nsd + 203);
    const auto *nsd_204 = buffer.data(nsd + 204);
    const auto *nsd_207 = buffer.data(nsd + 207);
    const auto *nsd_209 = buffer.data(nsd + 209);
    const auto *nsd_210 = buffer.data(nsd + 210);
    const auto *nsd_213 = buffer.data(nsd + 213);
    const auto *nsd_215 = buffer.data(nsd + 215);
    const auto *nsd_216 = buffer.data(nsd + 216);
    const auto *nsd_219 = buffer.data(nsd + 219);
    const auto *nsd_221 = buffer.data(nsd + 221);
    const auto *nsd_222 = buffer.data(nsd + 222);
    const auto *nsd_225 = buffer.data(nsd + 225);
    const auto *nsd_227 = buffer.data(nsd + 227);
    const auto *nsd_228 = buffer.data(nsd + 228);
    const auto *nsd_231 = buffer.data(nsd + 231);
    const auto *nsd_232 = buffer.data(nsd + 232);
    const auto *nsd_233 = buffer.data(nsd + 233);
    const auto *nsd_234 = buffer.data(nsd + 234);
    const auto *nsd_237 = buffer.data(nsd + 237);
    const auto *nsd_238 = buffer.data(nsd + 238);
    const auto *nsd_239 = buffer.data(nsd + 239);
    const auto *nsd_240 = buffer.data(nsd + 240);
    const auto *nsd_243 = buffer.data(nsd + 243);
    const auto *nsd_244 = buffer.data(nsd + 244);
    const auto *nsd_245 = buffer.data(nsd + 245);
    const auto *nsd_246 = buffer.data(nsd + 246);
    const auto *nsd_249 = buffer.data(nsd + 249);
    const auto *nsd_250 = buffer.data(nsd + 250);
    const auto *nsd_251 = buffer.data(nsd + 251);
    const auto *nsd_252 = buffer.data(nsd + 252);
    const auto *nsd_255 = buffer.data(nsd + 255);
    const auto *nsd_256 = buffer.data(nsd + 256);
    const auto *nsd_257 = buffer.data(nsd + 257);
    const auto *nsd_261 = buffer.data(nsd + 261);
    const auto *nsd_262 = buffer.data(nsd + 262);
    const auto *nsd_263 = buffer.data(nsd + 263);
    const auto *nsd_264 = buffer.data(nsd + 264);
    const auto *nsd_267 = buffer.data(nsd + 267);
    const auto *nsd_269 = buffer.data(nsd + 269);
    const auto *nsd_270 = buffer.data(nsd + 270);
    const auto *nsd_273 = buffer.data(nsd + 273);
    const auto *nsd_275 = buffer.data(nsd + 275);
    const auto *nsd_279 = buffer.data(nsd + 279);
    const auto *nsd_280 = buffer.data(nsd + 280);
    const auto *nsd_281 = buffer.data(nsd + 281);
    const auto *nsd_282 = buffer.data(nsd + 282);
    const auto *nsd_285 = buffer.data(nsd + 285);
    const auto *nsd_286 = buffer.data(nsd + 286);
    const auto *nsd_287 = buffer.data(nsd + 287);
    const auto *nsd_288 = buffer.data(nsd + 288);
    const auto *nsd_291 = buffer.data(nsd + 291);
    const auto *nsd_292 = buffer.data(nsd + 292);
    const auto *nsd_293 = buffer.data(nsd + 293);
    const auto *nsd_294 = buffer.data(nsd + 294);
    const auto *nsd_297 = buffer.data(nsd + 297);
    const auto *nsd_298 = buffer.data(nsd + 298);
    const auto *nsd_299 = buffer.data(nsd + 299);

    const auto *nsf1_350 = buffer.data(nsf1 + 350);
    const auto *nsf1_359 = buffer.data(nsf1 + 359);
    const auto *nsf1_360 = buffer.data(nsf1 + 360);
    const auto *nsf1_366 = buffer.data(nsf1 + 366);

    const auto *osp0_114 = buffer.data(osp0 + 114);
    const auto *osp0_115 = buffer.data(osp0 + 115);
    const auto *osp0_116 = buffer.data(osp0 + 116);
    const auto *osp0_117 = buffer.data(osp0 + 117);
    const auto *osp0_118 = buffer.data(osp0 + 118);
    const auto *osp0_119 = buffer.data(osp0 + 119);
    const auto *osp0_120 = buffer.data(osp0 + 120);
    const auto *osp0_121 = buffer.data(osp0 + 121);
    const auto *osp0_122 = buffer.data(osp0 + 122);
    const auto *osp0_123 = buffer.data(osp0 + 123);
    const auto *osp0_124 = buffer.data(osp0 + 124);
    const auto *osp0_125 = buffer.data(osp0 + 125);
    const auto *osp0_126 = buffer.data(osp0 + 126);
    const auto *osp0_127 = buffer.data(osp0 + 127);
    const auto *osp0_128 = buffer.data(osp0 + 128);
    const auto *osp0_130 = buffer.data(osp0 + 130);
    const auto *osp0_132 = buffer.data(osp0 + 132);
    const auto *osp0_133 = buffer.data(osp0 + 133);
    const auto *osp0_134 = buffer.data(osp0 + 134);
    const auto *osp0_135 = buffer.data(osp0 + 135);
    const auto *osp0_136 = buffer.data(osp0 + 136);
    const auto *osp0_137 = buffer.data(osp0 + 137);
    const auto *osp0_140 = buffer.data(osp0 + 140);
    const auto *osp0_141 = buffer.data(osp0 + 141);
    const auto *osp0_142 = buffer.data(osp0 + 142);
    const auto *osp0_143 = buffer.data(osp0 + 143);
    const auto *osp0_144 = buffer.data(osp0 + 144);
    const auto *osp0_145 = buffer.data(osp0 + 145);
    const auto *osp0_146 = buffer.data(osp0 + 146);
    const auto *osp0_147 = buffer.data(osp0 + 147);
    const auto *osp0_148 = buffer.data(osp0 + 148);
    const auto *osp0_149 = buffer.data(osp0 + 149);

    const auto *osp1_114 = buffer.data(osp1 + 114);
    const auto *osp1_115 = buffer.data(osp1 + 115);
    const auto *osp1_116 = buffer.data(osp1 + 116);
    const auto *osp1_117 = buffer.data(osp1 + 117);
    const auto *osp1_118 = buffer.data(osp1 + 118);
    const auto *osp1_119 = buffer.data(osp1 + 119);
    const auto *osp1_120 = buffer.data(osp1 + 120);
    const auto *osp1_121 = buffer.data(osp1 + 121);
    const auto *osp1_122 = buffer.data(osp1 + 122);
    const auto *osp1_123 = buffer.data(osp1 + 123);
    const auto *osp1_124 = buffer.data(osp1 + 124);
    const auto *osp1_125 = buffer.data(osp1 + 125);
    const auto *osp1_126 = buffer.data(osp1 + 126);
    const auto *osp1_127 = buffer.data(osp1 + 127);
    const auto *osp1_128 = buffer.data(osp1 + 128);
    const auto *osp1_130 = buffer.data(osp1 + 130);
    const auto *osp1_132 = buffer.data(osp1 + 132);
    const auto *osp1_133 = buffer.data(osp1 + 133);
    const auto *osp1_134 = buffer.data(osp1 + 134);
    const auto *osp1_135 = buffer.data(osp1 + 135);
    const auto *osp1_136 = buffer.data(osp1 + 136);
    const auto *osp1_137 = buffer.data(osp1 + 137);
    const auto *osp1_140 = buffer.data(osp1 + 140);
    const auto *osp1_141 = buffer.data(osp1 + 141);
    const auto *osp1_142 = buffer.data(osp1 + 142);
    const auto *osp1_143 = buffer.data(osp1 + 143);
    const auto *osp1_144 = buffer.data(osp1 + 144);
    const auto *osp1_145 = buffer.data(osp1 + 145);
    const auto *osp1_146 = buffer.data(osp1 + 146);
    const auto *osp1_147 = buffer.data(osp1 + 147);
    const auto *osp1_148 = buffer.data(osp1 + 148);
    const auto *osp1_149 = buffer.data(osp1 + 149);

    const auto *osd_228 = buffer.data(osd + 228);
    const auto *osd_231 = buffer.data(osd + 231);
    const auto *osd_232 = buffer.data(osd + 232);
    const auto *osd_233 = buffer.data(osd + 233);
    const auto *osd_234 = buffer.data(osd + 234);
    const auto *osd_237 = buffer.data(osd + 237);
    const auto *osd_238 = buffer.data(osd + 238);
    const auto *osd_239 = buffer.data(osd + 239);
    const auto *osd_240 = buffer.data(osd + 240);
    const auto *osd_243 = buffer.data(osd + 243);
    const auto *osd_244 = buffer.data(osd + 244);
    const auto *osd_245 = buffer.data(osd + 245);
    const auto *osd_246 = buffer.data(osd + 246);
    const auto *osd_249 = buffer.data(osd + 249);
    const auto *osd_250 = buffer.data(osd + 250);
    const auto *osd_251 = buffer.data(osd + 251);
    const auto *osd_252 = buffer.data(osd + 252);
    const auto *osd_255 = buffer.data(osd + 255);
    const auto *osd_256 = buffer.data(osd + 256);
    const auto *osd_257 = buffer.data(osd + 257);
    const auto *osd_258 = buffer.data(osd + 258);
    const auto *osd_261 = buffer.data(osd + 261);
    const auto *osd_262 = buffer.data(osd + 262);
    const auto *osd_263 = buffer.data(osd + 263);
    const auto *osd_264 = buffer.data(osd + 264);
    const auto *osd_266 = buffer.data(osd + 266);
    const auto *osd_267 = buffer.data(osd + 267);
    const auto *osd_268 = buffer.data(osd + 268);
    const auto *osd_269 = buffer.data(osd + 269);
    const auto *osd_270 = buffer.data(osd + 270);
    const auto *osd_271 = buffer.data(osd + 271);
    const auto *osd_273 = buffer.data(osd + 273);
    const auto *osd_275 = buffer.data(osd + 275);
    const auto *osd_276 = buffer.data(osd + 276);
    const auto *osd_279 = buffer.data(osd + 279);
    const auto *osd_280 = buffer.data(osd + 280);
    const auto *osd_281 = buffer.data(osd + 281);
    const auto *osd_282 = buffer.data(osd + 282);
    const auto *osd_285 = buffer.data(osd + 285);
    const auto *osd_286 = buffer.data(osd + 286);
    const auto *osd_287 = buffer.data(osd + 287);
    const auto *osd_288 = buffer.data(osd + 288);
    const auto *osd_291 = buffer.data(osd + 291);
    const auto *osd_292 = buffer.data(osd + 292);
    const auto *osd_293 = buffer.data(osd + 293);
    const auto *osd_294 = buffer.data(osd + 294);
    const auto *osd_297 = buffer.data(osd + 297);
    const auto *osd_298 = buffer.data(osd + 298);
    const auto *osd_299 = buffer.data(osd + 299);

#pragma omp simd aligned(t_380, t_381, t_382, t_383, pc_x, pc_y, pc_z, nsd_174, nsd_180, \
                         nsd_228, nsd_231, osp0_114, osp1_114, osd_228, \
                         osd_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = f_12 * nsd_228[k]
                   + f_1 * osp0_114[k]
                   - f_2 * osp1_114[k]
                   + f_3 * pc_x[k] * osd_228[k];

        t_381[k] = f_15 * nsd_180[k]
                   + f_3 * pc_y[k] * osd_228[k];

        t_382[k] = f_10 * nsd_174[k]
                   + f_3 * pc_z[k] * osd_228[k];

        t_383[k] = f_12 * nsd_231[k]
                   + f_3 * pc_x[k] * osd_231[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, t_387, pc_x, pc_y, pc_z, nsd_177, nsd_183, \
                         nsd_232, nsd_233, osp0_115, osp1_115, osd_231, osd_232, \
                         osd_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = f_12 * nsd_232[k]
                   + f_3 * pc_x[k] * osd_232[k];

        t_385[k] = f_12 * nsd_233[k]
                   + f_3 * pc_x[k] * osd_233[k];

        t_386[k] = f_15 * nsd_183[k]
                   + f_1 * osp0_115[k]
                   - f_2 * osp1_115[k]
                   + f_3 * pc_y[k] * osd_231[k];

        t_387[k] = f_10 * nsd_177[k]
                   + f_3 * pc_z[k] * osd_231[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, pc_x, pc_y, pc_z, nsd_179, nsd_185, nsd_234, \
                         osp0_116, osp0_117, osp1_116, osp1_117, osd_233, \
                         osd_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_15 * nsd_185[k]
                   + f_3 * pc_y[k] * osd_233[k];

        t_389[k] = f_10 * nsd_179[k]
                   + f_1 * osp0_116[k]
                   - f_2 * osp1_116[k]
                   + f_3 * pc_z[k] * osd_233[k];

        t_390[k] = f_12 * nsd_234[k]
                   + f_1 * osp0_117[k]
                   - f_2 * osp1_117[k]
                   + f_3 * pc_x[k] * osd_234[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, pc_x, pc_y, pc_z, nsd_180, nsd_186, \
                         nsd_237, nsd_238, osd_234, osd_237, osd_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_16 * nsd_186[k]
                   + f_3 * pc_y[k] * osd_234[k];

        t_392[k] = f_12 * nsd_180[k]
                   + f_3 * pc_z[k] * osd_234[k];

        t_393[k] = f_12 * nsd_237[k]
                   + f_3 * pc_x[k] * osd_237[k];

        t_394[k] = f_12 * nsd_238[k]
                   + f_3 * pc_x[k] * osd_238[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, pc_x, pc_y, pc_z, nsd_183, nsd_189, \
                         nsd_191, nsd_239, osp0_118, osp1_118, osd_237, \
                         osd_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = f_12 * nsd_239[k]
                   + f_3 * pc_x[k] * osd_239[k];

        t_396[k] = f_16 * nsd_189[k]
                   + f_1 * osp0_118[k]
                   - f_2 * osp1_118[k]
                   + f_3 * pc_y[k] * osd_237[k];

        t_397[k] = f_12 * nsd_183[k]
                   + f_3 * pc_z[k] * osd_237[k];

        t_398[k] = f_16 * nsd_191[k]
                   + f_3 * pc_y[k] * osd_239[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, pc_x, pc_y, pc_z, nsd_185, nsd_192, nsd_240, \
                         osp0_119, osp0_120, osp1_119, osp1_120, osd_239, \
                         osd_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_399[k] = f_12 * nsd_185[k]
                   + f_1 * osp0_119[k]
                   - f_2 * osp1_119[k]
                   + f_3 * pc_z[k] * osd_239[k];

        t_400[k] = f_12 * nsd_240[k]
                   + f_1 * osp0_120[k]
                   - f_2 * osp1_120[k]
                   + f_3 * pc_x[k] * osd_240[k];

        t_401[k] = f_14 * nsd_192[k]
                   + f_3 * pc_y[k] * osd_240[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, t_405, pc_x, pc_z, nsd_186, nsd_243, nsd_244, \
                         nsd_245, osd_240, osd_243, osd_244, osd_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = f_14 * nsd_186[k]
                   + f_3 * pc_z[k] * osd_240[k];

        t_403[k] = f_12 * nsd_243[k]
                   + f_3 * pc_x[k] * osd_243[k];

        t_404[k] = f_12 * nsd_244[k]
                   + f_3 * pc_x[k] * osd_244[k];

        t_405[k] = f_12 * nsd_245[k]
                   + f_3 * pc_x[k] * osd_245[k];
    }

#pragma omp simd aligned(t_406, t_407, t_408, t_409, pc_y, pc_z, nsd_189, nsd_191, nsd_195, \
                         nsd_197, osp0_121, osp0_122, osp1_121, osp1_122, osd_243, \
                         osd_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_406[k] = f_14 * nsd_195[k]
                   + f_1 * osp0_121[k]
                   - f_2 * osp1_121[k]
                   + f_3 * pc_y[k] * osd_243[k];

        t_407[k] = f_14 * nsd_189[k]
                   + f_3 * pc_z[k] * osd_243[k];

        t_408[k] = f_14 * nsd_197[k]
                   + f_3 * pc_y[k] * osd_245[k];

        t_409[k] = f_14 * nsd_191[k]
                   + f_1 * osp0_122[k]
                   - f_2 * osp1_122[k]
                   + f_3 * pc_z[k] * osd_245[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, pc_x, pc_y, pc_z, nsd_192, nsd_198, \
                         nsd_246, nsd_249, osp0_123, osp1_123, osd_246, \
                         osd_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_12 * nsd_246[k]
                   + f_1 * osp0_123[k]
                   - f_2 * osp1_123[k]
                   + f_3 * pc_x[k] * osd_246[k];

        t_411[k] = f_12 * nsd_198[k]
                   + f_3 * pc_y[k] * osd_246[k];

        t_412[k] = f_16 * nsd_192[k]
                   + f_3 * pc_z[k] * osd_246[k];

        t_413[k] = f_12 * nsd_249[k]
                   + f_3 * pc_x[k] * osd_249[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, pc_x, pc_y, pc_z, nsd_195, nsd_201, \
                         nsd_250, nsd_251, osp0_124, osp1_124, osd_249, osd_250, \
                         osd_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_12 * nsd_250[k]
                   + f_3 * pc_x[k] * osd_250[k];

        t_415[k] = f_12 * nsd_251[k]
                   + f_3 * pc_x[k] * osd_251[k];

        t_416[k] = f_12 * nsd_201[k]
                   + f_1 * osp0_124[k]
                   - f_2 * osp1_124[k]
                   + f_3 * pc_y[k] * osd_249[k];

        t_417[k] = f_16 * nsd_195[k]
                   + f_3 * pc_z[k] * osd_249[k];
    }

#pragma omp simd aligned(t_418, t_419, t_420, pc_x, pc_y, pc_z, nsd_197, nsd_203, nsd_252, \
                         osp0_125, osp0_126, osp1_125, osp1_126, osd_251, \
                         osd_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_418[k] = f_12 * nsd_203[k]
                   + f_3 * pc_y[k] * osd_251[k];

        t_419[k] = f_16 * nsd_197[k]
                   + f_1 * osp0_125[k]
                   - f_2 * osp1_125[k]
                   + f_3 * pc_z[k] * osd_251[k];

        t_420[k] = f_12 * nsd_252[k]
                   + f_1 * osp0_126[k]
                   - f_2 * osp1_126[k]
                   + f_3 * pc_x[k] * osd_252[k];
    }

#pragma omp simd aligned(t_421, t_422, t_423, t_424, pc_x, pc_y, pc_z, nsd_198, nsd_204, \
                         nsd_255, nsd_256, osd_252, osd_255, osd_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_421[k] = f_10 * nsd_204[k]
                   + f_3 * pc_y[k] * osd_252[k];

        t_422[k] = f_15 * nsd_198[k]
                   + f_3 * pc_z[k] * osd_252[k];

        t_423[k] = f_12 * nsd_255[k]
                   + f_3 * pc_x[k] * osd_255[k];

        t_424[k] = f_12 * nsd_256[k]
                   + f_3 * pc_x[k] * osd_256[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, pc_x, pc_y, pc_z, nsd_201, nsd_207, \
                         nsd_209, nsd_257, osp0_127, osp1_127, osd_255, \
                         osd_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = f_12 * nsd_257[k]
                   + f_3 * pc_x[k] * osd_257[k];

        t_426[k] = f_10 * nsd_207[k]
                   + f_1 * osp0_127[k]
                   - f_2 * osp1_127[k]
                   + f_3 * pc_y[k] * osd_255[k];

        t_427[k] = f_15 * nsd_201[k]
                   + f_3 * pc_z[k] * osd_255[k];

        t_428[k] = f_10 * nsd_209[k]
                   + f_3 * pc_y[k] * osd_257[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, pa_y, pc_y, pc_z, nsf0_350, nsd_203, \
                         nsd_204, nsd_210, nsf1_350, osp0_128, osp1_128, osd_257, \
                         osd_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_15 * nsd_203[k]
                   + f_1 * osp0_128[k]
                   - f_2 * osp1_128[k]
                   + f_3 * pc_z[k] * osd_257[k];

        t_430[k] = pa_y[k] * nsf0_350[k]
                   - f_4 * pc_y[k] * nsf1_350[k];

        t_431[k] = f_5 * nsd_210[k]
                   + f_3 * pc_y[k] * osd_258[k];

        t_432[k] = f_13 * nsd_204[k]
                   + f_3 * pc_z[k] * osd_258[k];
    }

#pragma omp simd aligned(t_433, t_434, t_435, t_436, pc_x, pc_y, nsd_213, nsd_261, nsd_262, \
                         nsd_263, osp0_130, osp1_130, osd_261, osd_262, \
                         osd_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_433[k] = f_12 * nsd_261[k]
                   + f_3 * pc_x[k] * osd_261[k];

        t_434[k] = f_12 * nsd_262[k]
                   + f_3 * pc_x[k] * osd_262[k];

        t_435[k] = f_12 * nsd_263[k]
                   + f_3 * pc_x[k] * osd_263[k];

        t_436[k] = f_5 * nsd_213[k]
                   + f_1 * osp0_130[k]
                   - f_2 * osp1_130[k]
                   + f_3 * pc_y[k] * osd_261[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, pa_y, pc_y, pc_z, nsf0_359, nsd_207, nsd_215, \
                         nsf1_359, osd_261, osd_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_13 * nsd_207[k]
                   + f_3 * pc_z[k] * osd_261[k];

        t_438[k] = f_5 * nsd_215[k]
                   + f_3 * pc_y[k] * osd_263[k];

        t_439[k] = pa_y[k] * nsf0_359[k]
                   - f_4 * pc_y[k] * nsf1_359[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, pc_x, pc_y, pc_z, nsd_210, \
                         nsd_264, nsd_267, osp0_132, osp1_132, osd_264, osd_266, \
                         osd_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_12 * nsd_264[k]
                   + f_1 * osp0_132[k]
                   - f_2 * osp1_132[k]
                   + f_3 * pc_x[k] * osd_264[k];

        t_441[k] = f_3 * pc_y[k] * osd_264[k];

        t_442[k] = f_11 * nsd_210[k]
                   + f_3 * pc_z[k] * osd_264[k];

        t_443[k] = f_12 * nsd_267[k]
                   + f_3 * pc_x[k] * osd_267[k];

        t_444[k] = f_3 * pc_y[k] * osd_266[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, pc_x, pc_y, nsd_269, osp0_133, osp0_134, \
                         osp1_133, osp1_134, osd_267, osd_268, \
                         osd_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_12 * nsd_269[k]
                   + f_3 * pc_x[k] * osd_269[k];

        t_446[k] = f_1 * osp0_133[k]
                   - f_2 * osp1_133[k]
                   + f_3 * pc_y[k] * osd_267[k];

        t_447[k] = f_7 * osp0_134[k]
                   - f_8 * osp1_134[k]
                   + f_3 * pc_y[k] * osd_268[k];

        t_448[k] = f_3 * pc_y[k] * osd_269[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, pc_x, pc_y, pc_z, nsd_215, nsd_216, \
                         nsd_270, osp0_134, osp0_135, osp1_134, osp1_135, osd_269, \
                         osd_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_11 * nsd_215[k]
                   + f_1 * osp0_134[k]
                   - f_2 * osp1_134[k]
                   + f_3 * pc_z[k] * osd_269[k];

        t_450[k] = f_10 * nsd_270[k]
                   + f_1 * osp0_135[k]
                   - f_2 * osp1_135[k]
                   + f_3 * pc_x[k] * osd_270[k];

        t_451[k] = f_9 * nsd_216[k]
                   + f_3 * pc_y[k] * osd_270[k];

        t_452[k] = f_3 * pc_z[k] * osd_270[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, t_456, t_457, pc_x, pc_y, pc_z, nsd_219, \
                         nsd_273, nsd_275, osp0_136, osp1_136, osd_271, osd_273, \
                         osd_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_10 * nsd_273[k]
                   + f_3 * pc_x[k] * osd_273[k];

        t_454[k] = f_3 * pc_z[k] * osd_271[k];

        t_455[k] = f_10 * nsd_275[k]
                   + f_3 * pc_x[k] * osd_275[k];

        t_456[k] = f_9 * nsd_219[k]
                   + f_1 * osp0_136[k]
                   - f_2 * osp1_136[k]
                   + f_3 * pc_y[k] * osd_273[k];

        t_457[k] = f_3 * pc_z[k] * osd_273[k];
    }

#pragma omp simd aligned(t_458, t_459, t_460, t_461, pa_z, pc_y, pc_z, nsf0_360, nsd_221, \
                         nsd_222, nsf1_360, osp0_137, osp1_137, osd_275, \
                         osd_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_458[k] = f_9 * nsd_221[k]
                   + f_3 * pc_y[k] * osd_275[k];

        t_459[k] = f_1 * osp0_137[k]
                   - f_2 * osp1_137[k]
                   + f_3 * pc_z[k] * osd_275[k];

        t_460[k] = pa_z[k] * nsf0_360[k]
                   - f_4 * pc_z[k] * nsf1_360[k];

        t_461[k] = f_11 * nsd_222[k]
                   + f_3 * pc_y[k] * osd_276[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, t_465, pc_x, pc_z, nsd_216, nsd_279, nsd_280, \
                         nsd_281, osd_276, osd_279, osd_280, osd_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = f_5 * nsd_216[k]
                   + f_3 * pc_z[k] * osd_276[k];

        t_463[k] = f_10 * nsd_279[k]
                   + f_3 * pc_x[k] * osd_279[k];

        t_464[k] = f_10 * nsd_280[k]
                   + f_3 * pc_x[k] * osd_280[k];

        t_465[k] = f_10 * nsd_281[k]
                   + f_3 * pc_x[k] * osd_281[k];
    }

#pragma omp simd aligned(t_466, t_467, t_468, t_469, pa_z, pc_y, pc_z, nsf0_366, nsd_219, \
                         nsd_221, nsd_227, nsf1_366, osp0_140, osp1_140, osd_279, \
                         osd_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_466[k] = pa_z[k] * nsf0_366[k]
                   - f_4 * pc_z[k] * nsf1_366[k];

        t_467[k] = f_5 * nsd_219[k]
                   + f_3 * pc_z[k] * osd_279[k];

        t_468[k] = f_11 * nsd_227[k]
                   + f_3 * pc_y[k] * osd_281[k];

        t_469[k] = f_5 * nsd_221[k]
                   + f_1 * osp0_140[k]
                   - f_2 * osp1_140[k]
                   + f_3 * pc_z[k] * osd_281[k];
    }

#pragma omp simd aligned(t_470, t_471, t_472, t_473, pc_x, pc_y, pc_z, nsd_222, nsd_228, \
                         nsd_282, nsd_285, osp0_141, osp1_141, osd_282, \
                         osd_285 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_470[k] = f_10 * nsd_282[k]
                   + f_1 * osp0_141[k]
                   - f_2 * osp1_141[k]
                   + f_3 * pc_x[k] * osd_282[k];

        t_471[k] = f_13 * nsd_228[k]
                   + f_3 * pc_y[k] * osd_282[k];

        t_472[k] = f_10 * nsd_222[k]
                   + f_3 * pc_z[k] * osd_282[k];

        t_473[k] = f_10 * nsd_285[k]
                   + f_3 * pc_x[k] * osd_285[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, pc_x, pc_y, pc_z, nsd_225, nsd_231, \
                         nsd_286, nsd_287, osp0_142, osp1_142, osd_285, osd_286, \
                         osd_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_10 * nsd_286[k]
                   + f_3 * pc_x[k] * osd_286[k];

        t_475[k] = f_10 * nsd_287[k]
                   + f_3 * pc_x[k] * osd_287[k];

        t_476[k] = f_13 * nsd_231[k]
                   + f_1 * osp0_142[k]
                   - f_2 * osp1_142[k]
                   + f_3 * pc_y[k] * osd_285[k];

        t_477[k] = f_10 * nsd_225[k]
                   + f_3 * pc_z[k] * osd_285[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, pc_x, pc_y, pc_z, nsd_227, nsd_233, nsd_288, \
                         osp0_143, osp0_144, osp1_143, osp1_144, osd_287, \
                         osd_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = f_13 * nsd_233[k]
                   + f_3 * pc_y[k] * osd_287[k];

        t_479[k] = f_10 * nsd_227[k]
                   + f_1 * osp0_143[k]
                   - f_2 * osp1_143[k]
                   + f_3 * pc_z[k] * osd_287[k];

        t_480[k] = f_10 * nsd_288[k]
                   + f_1 * osp0_144[k]
                   - f_2 * osp1_144[k]
                   + f_3 * pc_x[k] * osd_288[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, t_484, pc_x, pc_y, pc_z, nsd_228, nsd_234, \
                         nsd_291, nsd_292, osd_288, osd_291, osd_292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = f_15 * nsd_234[k]
                   + f_3 * pc_y[k] * osd_288[k];

        t_482[k] = f_12 * nsd_228[k]
                   + f_3 * pc_z[k] * osd_288[k];

        t_483[k] = f_10 * nsd_291[k]
                   + f_3 * pc_x[k] * osd_291[k];

        t_484[k] = f_10 * nsd_292[k]
                   + f_3 * pc_x[k] * osd_292[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, t_488, pc_x, pc_y, pc_z, nsd_231, nsd_237, \
                         nsd_239, nsd_293, osp0_145, osp1_145, osd_291, \
                         osd_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_10 * nsd_293[k]
                   + f_3 * pc_x[k] * osd_293[k];

        t_486[k] = f_15 * nsd_237[k]
                   + f_1 * osp0_145[k]
                   - f_2 * osp1_145[k]
                   + f_3 * pc_y[k] * osd_291[k];

        t_487[k] = f_12 * nsd_231[k]
                   + f_3 * pc_z[k] * osd_291[k];

        t_488[k] = f_15 * nsd_239[k]
                   + f_3 * pc_y[k] * osd_293[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, pc_x, pc_y, pc_z, nsd_233, nsd_240, nsd_294, \
                         osp0_146, osp0_147, osp1_146, osp1_147, osd_293, \
                         osd_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = f_12 * nsd_233[k]
                   + f_1 * osp0_146[k]
                   - f_2 * osp1_146[k]
                   + f_3 * pc_z[k] * osd_293[k];

        t_490[k] = f_10 * nsd_294[k]
                   + f_1 * osp0_147[k]
                   - f_2 * osp1_147[k]
                   + f_3 * pc_x[k] * osd_294[k];

        t_491[k] = f_16 * nsd_240[k]
                   + f_3 * pc_y[k] * osd_294[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, pc_x, pc_z, nsd_234, nsd_297, nsd_298, \
                         nsd_299, osd_294, osd_297, osd_298, osd_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = f_14 * nsd_234[k]
                   + f_3 * pc_z[k] * osd_294[k];

        t_493[k] = f_10 * nsd_297[k]
                   + f_3 * pc_x[k] * osd_297[k];

        t_494[k] = f_10 * nsd_298[k]
                   + f_3 * pc_x[k] * osd_298[k];

        t_495[k] = f_10 * nsd_299[k]
                   + f_3 * pc_x[k] * osd_299[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, t_499, pc_y, pc_z, nsd_237, nsd_239, nsd_243, \
                         nsd_245, osp0_148, osp0_149, osp1_148, osp1_149, osd_297, \
                         osd_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_16 * nsd_243[k]
                   + f_1 * osp0_148[k]
                   - f_2 * osp1_148[k]
                   + f_3 * pc_y[k] * osd_297[k];

        t_497[k] = f_14 * nsd_237[k]
                   + f_3 * pc_z[k] * osd_297[k];

        t_498[k] = f_16 * nsd_245[k]
                   + f_3 * pc_y[k] * osd_299[k];

        t_499[k] = f_14 * nsd_239[k]
                   + f_1 * osp0_149[k]
                   - f_2 * osp1_149[k]
                   + f_3 * pc_z[k] * osd_299[k];
    }
}

static auto
compute_prim_osf_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsf0,
                                                          const size_t nsd, const size_t nsf1,
                                                          const size_t osp0, const size_t osp1,
                                                          const size_t osd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 5.0 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_9 = 4.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 4.0 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 3.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsf0_440 = buffer.data(nsf0 + 440);
    const auto *nsf0_449 = buffer.data(nsf0 + 449);
    const auto *nsf0_450 = buffer.data(nsf0 + 450);
    const auto *nsf0_550 = buffer.data(nsf0 + 550);
    const auto *nsf0_556 = buffer.data(nsf0 + 556);
    const auto *nsf0_559 = buffer.data(nsf0 + 559);
    const auto *nsf0_566 = buffer.data(nsf0 + 566);
    const auto *nsf0_569 = buffer.data(nsf0 + 569);
    const auto *nsf0_570 = buffer.data(nsf0 + 570);
    const auto *nsf0_576 = buffer.data(nsf0 + 576);
    const auto *nsf0_579 = buffer.data(nsf0 + 579);
    const auto *nsf0_580 = buffer.data(nsf0 + 580);
    const auto *nsf0_586 = buffer.data(nsf0 + 586);
    const auto *nsf0_589 = buffer.data(nsf0 + 589);
    const auto *nsf0_590 = buffer.data(nsf0 + 590);
    const auto *nsf0_596 = buffer.data(nsf0 + 596);
    const auto *nsf0_599 = buffer.data(nsf0 + 599);
    const auto *nsf0_600 = buffer.data(nsf0 + 600);
    const auto *nsf0_606 = buffer.data(nsf0 + 606);
    const auto *nsf0_609 = buffer.data(nsf0 + 609);
    const auto *nsf0_610 = buffer.data(nsf0 + 610);
    const auto *nsf0_616 = buffer.data(nsf0 + 616);
    const auto *nsf0_619 = buffer.data(nsf0 + 619);
    const auto *nsf0_620 = buffer.data(nsf0 + 620);

    const auto *nsd_240 = buffer.data(nsd + 240);
    const auto *nsd_243 = buffer.data(nsd + 243);
    const auto *nsd_245 = buffer.data(nsd + 245);
    const auto *nsd_246 = buffer.data(nsd + 246);
    const auto *nsd_249 = buffer.data(nsd + 249);
    const auto *nsd_251 = buffer.data(nsd + 251);
    const auto *nsd_252 = buffer.data(nsd + 252);
    const auto *nsd_255 = buffer.data(nsd + 255);
    const auto *nsd_257 = buffer.data(nsd + 257);
    const auto *nsd_258 = buffer.data(nsd + 258);
    const auto *nsd_261 = buffer.data(nsd + 261);
    const auto *nsd_263 = buffer.data(nsd + 263);
    const auto *nsd_264 = buffer.data(nsd + 264);
    const auto *nsd_267 = buffer.data(nsd + 267);
    const auto *nsd_269 = buffer.data(nsd + 269);
    const auto *nsd_270 = buffer.data(nsd + 270);
    const auto *nsd_273 = buffer.data(nsd + 273);
    const auto *nsd_275 = buffer.data(nsd + 275);
    const auto *nsd_276 = buffer.data(nsd + 276);
    const auto *nsd_279 = buffer.data(nsd + 279);
    const auto *nsd_281 = buffer.data(nsd + 281);
    const auto *nsd_282 = buffer.data(nsd + 282);
    const auto *nsd_285 = buffer.data(nsd + 285);
    const auto *nsd_287 = buffer.data(nsd + 287);
    const auto *nsd_288 = buffer.data(nsd + 288);
    const auto *nsd_291 = buffer.data(nsd + 291);
    const auto *nsd_293 = buffer.data(nsd + 293);
    const auto *nsd_294 = buffer.data(nsd + 294);
    const auto *nsd_297 = buffer.data(nsd + 297);
    const auto *nsd_299 = buffer.data(nsd + 299);
    const auto *nsd_300 = buffer.data(nsd + 300);
    const auto *nsd_303 = buffer.data(nsd + 303);
    const auto *nsd_304 = buffer.data(nsd + 304);
    const auto *nsd_305 = buffer.data(nsd + 305);
    const auto *nsd_306 = buffer.data(nsd + 306);
    const auto *nsd_309 = buffer.data(nsd + 309);
    const auto *nsd_310 = buffer.data(nsd + 310);
    const auto *nsd_311 = buffer.data(nsd + 311);
    const auto *nsd_312 = buffer.data(nsd + 312);
    const auto *nsd_315 = buffer.data(nsd + 315);
    const auto *nsd_316 = buffer.data(nsd + 316);
    const auto *nsd_317 = buffer.data(nsd + 317);
    const auto *nsd_321 = buffer.data(nsd + 321);
    const auto *nsd_322 = buffer.data(nsd + 322);
    const auto *nsd_323 = buffer.data(nsd + 323);
    const auto *nsd_324 = buffer.data(nsd + 324);
    const auto *nsd_327 = buffer.data(nsd + 327);
    const auto *nsd_329 = buffer.data(nsd + 329);
    const auto *nsd_330 = buffer.data(nsd + 330);
    const auto *nsd_333 = buffer.data(nsd + 333);
    const auto *nsd_335 = buffer.data(nsd + 335);
    const auto *nsd_339 = buffer.data(nsd + 339);
    const auto *nsd_340 = buffer.data(nsd + 340);
    const auto *nsd_341 = buffer.data(nsd + 341);
    const auto *nsd_342 = buffer.data(nsd + 342);
    const auto *nsd_345 = buffer.data(nsd + 345);
    const auto *nsd_346 = buffer.data(nsd + 346);
    const auto *nsd_347 = buffer.data(nsd + 347);
    const auto *nsd_348 = buffer.data(nsd + 348);
    const auto *nsd_351 = buffer.data(nsd + 351);
    const auto *nsd_352 = buffer.data(nsd + 352);
    const auto *nsd_353 = buffer.data(nsd + 353);
    const auto *nsd_354 = buffer.data(nsd + 354);
    const auto *nsd_357 = buffer.data(nsd + 357);
    const auto *nsd_358 = buffer.data(nsd + 358);
    const auto *nsd_359 = buffer.data(nsd + 359);
    const auto *nsd_360 = buffer.data(nsd + 360);
    const auto *nsd_363 = buffer.data(nsd + 363);
    const auto *nsd_364 = buffer.data(nsd + 364);
    const auto *nsd_365 = buffer.data(nsd + 365);
    const auto *nsd_366 = buffer.data(nsd + 366);
    const auto *nsd_369 = buffer.data(nsd + 369);
    const auto *nsd_370 = buffer.data(nsd + 370);
    const auto *nsd_371 = buffer.data(nsd + 371);
    const auto *nsd_372 = buffer.data(nsd + 372);
    const auto *nsd_375 = buffer.data(nsd + 375);
    const auto *nsd_376 = buffer.data(nsd + 376);
    const auto *nsd_377 = buffer.data(nsd + 377);

    const auto *nsf1_440 = buffer.data(nsf1 + 440);
    const auto *nsf1_449 = buffer.data(nsf1 + 449);
    const auto *nsf1_450 = buffer.data(nsf1 + 450);
    const auto *nsf1_550 = buffer.data(nsf1 + 550);
    const auto *nsf1_556 = buffer.data(nsf1 + 556);
    const auto *nsf1_559 = buffer.data(nsf1 + 559);
    const auto *nsf1_566 = buffer.data(nsf1 + 566);
    const auto *nsf1_569 = buffer.data(nsf1 + 569);
    const auto *nsf1_570 = buffer.data(nsf1 + 570);
    const auto *nsf1_576 = buffer.data(nsf1 + 576);
    const auto *nsf1_579 = buffer.data(nsf1 + 579);
    const auto *nsf1_580 = buffer.data(nsf1 + 580);
    const auto *nsf1_586 = buffer.data(nsf1 + 586);
    const auto *nsf1_589 = buffer.data(nsf1 + 589);
    const auto *nsf1_590 = buffer.data(nsf1 + 590);
    const auto *nsf1_596 = buffer.data(nsf1 + 596);
    const auto *nsf1_599 = buffer.data(nsf1 + 599);
    const auto *nsf1_600 = buffer.data(nsf1 + 600);
    const auto *nsf1_606 = buffer.data(nsf1 + 606);
    const auto *nsf1_609 = buffer.data(nsf1 + 609);
    const auto *nsf1_610 = buffer.data(nsf1 + 610);
    const auto *nsf1_616 = buffer.data(nsf1 + 616);
    const auto *nsf1_619 = buffer.data(nsf1 + 619);
    const auto *nsf1_620 = buffer.data(nsf1 + 620);

    const auto *osp0_150 = buffer.data(osp0 + 150);
    const auto *osp0_151 = buffer.data(osp0 + 151);
    const auto *osp0_152 = buffer.data(osp0 + 152);
    const auto *osp0_153 = buffer.data(osp0 + 153);
    const auto *osp0_154 = buffer.data(osp0 + 154);
    const auto *osp0_155 = buffer.data(osp0 + 155);
    const auto *osp0_156 = buffer.data(osp0 + 156);
    const auto *osp0_157 = buffer.data(osp0 + 157);
    const auto *osp0_158 = buffer.data(osp0 + 158);
    const auto *osp0_160 = buffer.data(osp0 + 160);
    const auto *osp0_162 = buffer.data(osp0 + 162);
    const auto *osp0_163 = buffer.data(osp0 + 163);
    const auto *osp0_164 = buffer.data(osp0 + 164);

    const auto *osp1_150 = buffer.data(osp1 + 150);
    const auto *osp1_151 = buffer.data(osp1 + 151);
    const auto *osp1_152 = buffer.data(osp1 + 152);
    const auto *osp1_153 = buffer.data(osp1 + 153);
    const auto *osp1_154 = buffer.data(osp1 + 154);
    const auto *osp1_155 = buffer.data(osp1 + 155);
    const auto *osp1_156 = buffer.data(osp1 + 156);
    const auto *osp1_157 = buffer.data(osp1 + 157);
    const auto *osp1_158 = buffer.data(osp1 + 158);
    const auto *osp1_160 = buffer.data(osp1 + 160);
    const auto *osp1_162 = buffer.data(osp1 + 162);
    const auto *osp1_163 = buffer.data(osp1 + 163);
    const auto *osp1_164 = buffer.data(osp1 + 164);

    const auto *osd_300 = buffer.data(osd + 300);
    const auto *osd_303 = buffer.data(osd + 303);
    const auto *osd_304 = buffer.data(osd + 304);
    const auto *osd_305 = buffer.data(osd + 305);
    const auto *osd_306 = buffer.data(osd + 306);
    const auto *osd_309 = buffer.data(osd + 309);
    const auto *osd_310 = buffer.data(osd + 310);
    const auto *osd_311 = buffer.data(osd + 311);
    const auto *osd_312 = buffer.data(osd + 312);
    const auto *osd_315 = buffer.data(osd + 315);
    const auto *osd_316 = buffer.data(osd + 316);
    const auto *osd_317 = buffer.data(osd + 317);
    const auto *osd_318 = buffer.data(osd + 318);
    const auto *osd_321 = buffer.data(osd + 321);
    const auto *osd_322 = buffer.data(osd + 322);
    const auto *osd_323 = buffer.data(osd + 323);
    const auto *osd_324 = buffer.data(osd + 324);
    const auto *osd_326 = buffer.data(osd + 326);
    const auto *osd_327 = buffer.data(osd + 327);
    const auto *osd_328 = buffer.data(osd + 328);
    const auto *osd_329 = buffer.data(osd + 329);
    const auto *osd_330 = buffer.data(osd + 330);
    const auto *osd_331 = buffer.data(osd + 331);
    const auto *osd_333 = buffer.data(osd + 333);
    const auto *osd_335 = buffer.data(osd + 335);
    const auto *osd_336 = buffer.data(osd + 336);
    const auto *osd_339 = buffer.data(osd + 339);
    const auto *osd_340 = buffer.data(osd + 340);
    const auto *osd_341 = buffer.data(osd + 341);
    const auto *osd_342 = buffer.data(osd + 342);
    const auto *osd_345 = buffer.data(osd + 345);
    const auto *osd_346 = buffer.data(osd + 346);
    const auto *osd_347 = buffer.data(osd + 347);
    const auto *osd_348 = buffer.data(osd + 348);
    const auto *osd_351 = buffer.data(osd + 351);
    const auto *osd_352 = buffer.data(osd + 352);
    const auto *osd_353 = buffer.data(osd + 353);
    const auto *osd_354 = buffer.data(osd + 354);
    const auto *osd_357 = buffer.data(osd + 357);
    const auto *osd_358 = buffer.data(osd + 358);
    const auto *osd_359 = buffer.data(osd + 359);
    const auto *osd_360 = buffer.data(osd + 360);
    const auto *osd_363 = buffer.data(osd + 363);
    const auto *osd_364 = buffer.data(osd + 364);
    const auto *osd_365 = buffer.data(osd + 365);
    const auto *osd_366 = buffer.data(osd + 366);
    const auto *osd_369 = buffer.data(osd + 369);
    const auto *osd_370 = buffer.data(osd + 370);
    const auto *osd_371 = buffer.data(osd + 371);
    const auto *osd_372 = buffer.data(osd + 372);
    const auto *osd_375 = buffer.data(osd + 375);
    const auto *osd_376 = buffer.data(osd + 376);
    const auto *osd_377 = buffer.data(osd + 377);

#pragma omp simd aligned(t_500, t_501, t_502, t_503, pc_x, pc_y, pc_z, nsd_240, nsd_246, \
                         nsd_300, nsd_303, osp0_150, osp1_150, osd_300, \
                         osd_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = f_10 * nsd_300[k]
                   + f_1 * osp0_150[k]
                   - f_2 * osp1_150[k]
                   + f_3 * pc_x[k] * osd_300[k];

        t_501[k] = f_14 * nsd_246[k]
                   + f_3 * pc_y[k] * osd_300[k];

        t_502[k] = f_16 * nsd_240[k]
                   + f_3 * pc_z[k] * osd_300[k];

        t_503[k] = f_10 * nsd_303[k]
                   + f_3 * pc_x[k] * osd_303[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, t_507, pc_x, pc_y, pc_z, nsd_243, nsd_249, \
                         nsd_304, nsd_305, osp0_151, osp1_151, osd_303, osd_304, \
                         osd_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_10 * nsd_304[k]
                   + f_3 * pc_x[k] * osd_304[k];

        t_505[k] = f_10 * nsd_305[k]
                   + f_3 * pc_x[k] * osd_305[k];

        t_506[k] = f_14 * nsd_249[k]
                   + f_1 * osp0_151[k]
                   - f_2 * osp1_151[k]
                   + f_3 * pc_y[k] * osd_303[k];

        t_507[k] = f_16 * nsd_243[k]
                   + f_3 * pc_z[k] * osd_303[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, pc_x, pc_y, pc_z, nsd_245, nsd_251, nsd_306, \
                         osp0_152, osp0_153, osp1_152, osp1_153, osd_305, \
                         osd_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = f_14 * nsd_251[k]
                   + f_3 * pc_y[k] * osd_305[k];

        t_509[k] = f_16 * nsd_245[k]
                   + f_1 * osp0_152[k]
                   - f_2 * osp1_152[k]
                   + f_3 * pc_z[k] * osd_305[k];

        t_510[k] = f_10 * nsd_306[k]
                   + f_1 * osp0_153[k]
                   - f_2 * osp1_153[k]
                   + f_3 * pc_x[k] * osd_306[k];
    }

#pragma omp simd aligned(t_511, t_512, t_513, t_514, pc_x, pc_y, pc_z, nsd_246, nsd_252, \
                         nsd_309, nsd_310, osd_306, osd_309, osd_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_511[k] = f_12 * nsd_252[k]
                   + f_3 * pc_y[k] * osd_306[k];

        t_512[k] = f_15 * nsd_246[k]
                   + f_3 * pc_z[k] * osd_306[k];

        t_513[k] = f_10 * nsd_309[k]
                   + f_3 * pc_x[k] * osd_309[k];

        t_514[k] = f_10 * nsd_310[k]
                   + f_3 * pc_x[k] * osd_310[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, pc_x, pc_y, pc_z, nsd_249, nsd_255, \
                         nsd_257, nsd_311, osp0_154, osp1_154, osd_309, \
                         osd_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = f_10 * nsd_311[k]
                   + f_3 * pc_x[k] * osd_311[k];

        t_516[k] = f_12 * nsd_255[k]
                   + f_1 * osp0_154[k]
                   - f_2 * osp1_154[k]
                   + f_3 * pc_y[k] * osd_309[k];

        t_517[k] = f_15 * nsd_249[k]
                   + f_3 * pc_z[k] * osd_309[k];

        t_518[k] = f_12 * nsd_257[k]
                   + f_3 * pc_y[k] * osd_311[k];
    }

#pragma omp simd aligned(t_519, t_520, t_521, pc_x, pc_y, pc_z, nsd_251, nsd_258, nsd_312, \
                         osp0_155, osp0_156, osp1_155, osp1_156, osd_311, \
                         osd_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_519[k] = f_15 * nsd_251[k]
                   + f_1 * osp0_155[k]
                   - f_2 * osp1_155[k]
                   + f_3 * pc_z[k] * osd_311[k];

        t_520[k] = f_10 * nsd_312[k]
                   + f_1 * osp0_156[k]
                   - f_2 * osp1_156[k]
                   + f_3 * pc_x[k] * osd_312[k];

        t_521[k] = f_10 * nsd_258[k]
                   + f_3 * pc_y[k] * osd_312[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, pc_x, pc_z, nsd_252, nsd_315, nsd_316, \
                         nsd_317, osd_312, osd_315, osd_316, osd_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_13 * nsd_252[k]
                   + f_3 * pc_z[k] * osd_312[k];

        t_523[k] = f_10 * nsd_315[k]
                   + f_3 * pc_x[k] * osd_315[k];

        t_524[k] = f_10 * nsd_316[k]
                   + f_3 * pc_x[k] * osd_316[k];

        t_525[k] = f_10 * nsd_317[k]
                   + f_3 * pc_x[k] * osd_317[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, t_529, pc_y, pc_z, nsd_255, nsd_257, nsd_261, \
                         nsd_263, osp0_157, osp0_158, osp1_157, osp1_158, osd_315, \
                         osd_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_10 * nsd_261[k]
                   + f_1 * osp0_157[k]
                   - f_2 * osp1_157[k]
                   + f_3 * pc_y[k] * osd_315[k];

        t_527[k] = f_13 * nsd_255[k]
                   + f_3 * pc_z[k] * osd_315[k];

        t_528[k] = f_10 * nsd_263[k]
                   + f_3 * pc_y[k] * osd_317[k];

        t_529[k] = f_13 * nsd_257[k]
                   + f_1 * osp0_158[k]
                   - f_2 * osp1_158[k]
                   + f_3 * pc_z[k] * osd_317[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, pa_y, pc_x, pc_y, pc_z, nsf0_440, \
                         nsd_258, nsd_264, nsd_321, nsf1_440, osd_318, \
                         osd_321 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = pa_y[k] * nsf0_440[k]
                   - f_4 * pc_y[k] * nsf1_440[k];

        t_531[k] = f_5 * nsd_264[k]
                   + f_3 * pc_y[k] * osd_318[k];

        t_532[k] = f_11 * nsd_258[k]
                   + f_3 * pc_z[k] * osd_318[k];

        t_533[k] = f_10 * nsd_321[k]
                   + f_3 * pc_x[k] * osd_321[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, t_537, pc_x, pc_y, pc_z, nsd_261, nsd_267, \
                         nsd_322, nsd_323, osp0_160, osp1_160, osd_321, osd_322, \
                         osd_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = f_10 * nsd_322[k]
                   + f_3 * pc_x[k] * osd_322[k];

        t_535[k] = f_10 * nsd_323[k]
                   + f_3 * pc_x[k] * osd_323[k];

        t_536[k] = f_5 * nsd_267[k]
                   + f_1 * osp0_160[k]
                   - f_2 * osp1_160[k]
                   + f_3 * pc_y[k] * osd_321[k];

        t_537[k] = f_11 * nsd_261[k]
                   + f_3 * pc_z[k] * osd_321[k];
    }

#pragma omp simd aligned(t_538, t_539, t_540, t_541, pa_y, pc_x, pc_y, nsf0_449, nsd_269, \
                         nsd_324, nsf1_449, osp0_162, osp1_162, osd_323, \
                         osd_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_538[k] = f_5 * nsd_269[k]
                   + f_3 * pc_y[k] * osd_323[k];

        t_539[k] = pa_y[k] * nsf0_449[k]
                   - f_4 * pc_y[k] * nsf1_449[k];

        t_540[k] = f_10 * nsd_324[k]
                   + f_1 * osp0_162[k]
                   - f_2 * osp1_162[k]
                   + f_3 * pc_x[k] * osd_324[k];

        t_541[k] = f_3 * pc_y[k] * osd_324[k];
    }

#pragma omp simd aligned(t_542, t_543, t_544, t_545, pc_x, pc_y, pc_z, nsd_264, nsd_327, \
                         nsd_329, osd_324, osd_326, osd_327, osd_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_542[k] = f_9 * nsd_264[k]
                   + f_3 * pc_z[k] * osd_324[k];

        t_543[k] = f_10 * nsd_327[k]
                   + f_3 * pc_x[k] * osd_327[k];

        t_544[k] = f_3 * pc_y[k] * osd_326[k];

        t_545[k] = f_10 * nsd_329[k]
                   + f_3 * pc_x[k] * osd_329[k];
    }

#pragma omp simd aligned(t_546, t_547, t_548, t_549, pc_y, pc_z, nsd_269, osp0_163, osp0_164, \
                         osp1_163, osp1_164, osd_327, osd_328, \
                         osd_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_546[k] = f_1 * osp0_163[k]
                   - f_2 * osp1_163[k]
                   + f_3 * pc_y[k] * osd_327[k];

        t_547[k] = f_7 * osp0_164[k]
                   - f_8 * osp1_164[k]
                   + f_3 * pc_y[k] * osd_328[k];

        t_548[k] = f_3 * pc_y[k] * osd_329[k];

        t_549[k] = f_9 * nsd_269[k]
                   + f_1 * osp0_164[k]
                   - f_2 * osp1_164[k]
                   + f_3 * pc_z[k] * osd_329[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, t_553, pa_x, pc_x, pc_y, pc_z, nsf0_550, \
                         nsd_270, nsd_330, nsd_333, nsf1_550, osd_330, \
                         osd_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = pa_x[k] * nsf0_550[k]
                   + f_12 * nsd_330[k]
                   - f_4 * pc_x[k] * nsf1_550[k];

        t_551[k] = f_6 * nsd_270[k]
                   + f_3 * pc_y[k] * osd_330[k];

        t_552[k] = f_3 * pc_z[k] * osd_330[k];

        t_553[k] = f_5 * nsd_333[k]
                   + f_3 * pc_x[k] * osd_333[k];
    }

#pragma omp simd aligned(t_554, t_555, t_556, t_557, t_558, pa_x, pc_x, pc_y, pc_z, nsf0_556, \
                         nsd_275, nsd_335, nsf1_556, osd_331, osd_333, \
                         osd_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_554[k] = f_3 * pc_z[k] * osd_331[k];

        t_555[k] = f_5 * nsd_335[k]
                   + f_3 * pc_x[k] * osd_335[k];

        t_556[k] = pa_x[k] * nsf0_556[k]
                   - f_4 * pc_x[k] * nsf1_556[k];

        t_557[k] = f_3 * pc_z[k] * osd_333[k];

        t_558[k] = f_6 * nsd_275[k]
                   + f_3 * pc_y[k] * osd_335[k];
    }

#pragma omp simd aligned(t_559, t_560, t_561, t_562, pa_x, pa_z, pc_x, pc_y, pc_z, nsf0_450, \
                         nsf0_559, nsd_270, nsd_276, nsf1_450, nsf1_559, \
                         osd_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_559[k] = pa_x[k] * nsf0_559[k]
                   - f_4 * pc_x[k] * nsf1_559[k];

        t_560[k] = pa_z[k] * nsf0_450[k]
                   - f_4 * pc_z[k] * nsf1_450[k];

        t_561[k] = f_9 * nsd_276[k]
                   + f_3 * pc_y[k] * osd_336[k];

        t_562[k] = f_5 * nsd_270[k]
                   + f_3 * pc_z[k] * osd_336[k];
    }

#pragma omp simd aligned(t_563, t_564, t_565, t_566, pa_x, pc_x, nsf0_566, nsd_339, nsd_340, \
                         nsd_341, nsf1_566, osd_339, osd_340, osd_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_563[k] = f_5 * nsd_339[k]
                   + f_3 * pc_x[k] * osd_339[k];

        t_564[k] = f_5 * nsd_340[k]
                   + f_3 * pc_x[k] * osd_340[k];

        t_565[k] = f_5 * nsd_341[k]
                   + f_3 * pc_x[k] * osd_341[k];

        t_566[k] = pa_x[k] * nsf0_566[k]
                   - f_4 * pc_x[k] * nsf1_566[k];
    }

#pragma omp simd aligned(t_567, t_568, t_569, pa_x, pc_x, pc_y, pc_z, nsf0_569, nsd_273, \
                         nsd_281, nsf1_569, osd_339, osd_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_567[k] = f_5 * nsd_273[k]
                   + f_3 * pc_z[k] * osd_339[k];

        t_568[k] = f_9 * nsd_281[k]
                   + f_3 * pc_y[k] * osd_341[k];

        t_569[k] = pa_x[k] * nsf0_569[k]
                   - f_4 * pc_x[k] * nsf1_569[k];
    }

#pragma omp simd aligned(t_570, t_571, t_572, t_573, pa_x, pc_x, pc_y, pc_z, nsf0_570, \
                         nsd_276, nsd_282, nsd_342, nsd_345, nsf1_570, osd_342, \
                         osd_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_570[k] = pa_x[k] * nsf0_570[k]
                   + f_12 * nsd_342[k]
                   - f_4 * pc_x[k] * nsf1_570[k];

        t_571[k] = f_11 * nsd_282[k]
                   + f_3 * pc_y[k] * osd_342[k];

        t_572[k] = f_10 * nsd_276[k]
                   + f_3 * pc_z[k] * osd_342[k];

        t_573[k] = f_5 * nsd_345[k]
                   + f_3 * pc_x[k] * osd_345[k];
    }

#pragma omp simd aligned(t_574, t_575, t_576, t_577, pa_x, pc_x, pc_z, nsf0_576, nsd_279, \
                         nsd_346, nsd_347, nsf1_576, osd_345, osd_346, \
                         osd_347 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_574[k] = f_5 * nsd_346[k]
                   + f_3 * pc_x[k] * osd_346[k];

        t_575[k] = f_5 * nsd_347[k]
                   + f_3 * pc_x[k] * osd_347[k];

        t_576[k] = pa_x[k] * nsf0_576[k]
                   - f_4 * pc_x[k] * nsf1_576[k];

        t_577[k] = f_10 * nsd_279[k]
                   + f_3 * pc_z[k] * osd_345[k];
    }

#pragma omp simd aligned(t_578, t_579, t_580, t_581, pa_x, pc_x, pc_y, nsf0_579, nsf0_580, \
                         nsd_287, nsd_288, nsd_348, nsf1_579, nsf1_580, osd_347, \
                         osd_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_578[k] = f_11 * nsd_287[k]
                   + f_3 * pc_y[k] * osd_347[k];

        t_579[k] = pa_x[k] * nsf0_579[k]
                   - f_4 * pc_x[k] * nsf1_579[k];

        t_580[k] = pa_x[k] * nsf0_580[k]
                   + f_12 * nsd_348[k]
                   - f_4 * pc_x[k] * nsf1_580[k];

        t_581[k] = f_13 * nsd_288[k]
                   + f_3 * pc_y[k] * osd_348[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, t_585, pc_x, pc_z, nsd_282, nsd_351, nsd_352, \
                         nsd_353, osd_348, osd_351, osd_352, osd_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = f_12 * nsd_282[k]
                   + f_3 * pc_z[k] * osd_348[k];

        t_583[k] = f_5 * nsd_351[k]
                   + f_3 * pc_x[k] * osd_351[k];

        t_584[k] = f_5 * nsd_352[k]
                   + f_3 * pc_x[k] * osd_352[k];

        t_585[k] = f_5 * nsd_353[k]
                   + f_3 * pc_x[k] * osd_353[k];
    }

#pragma omp simd aligned(t_586, t_587, t_588, t_589, pa_x, pc_x, pc_y, pc_z, nsf0_586, \
                         nsf0_589, nsd_285, nsd_293, nsf1_586, nsf1_589, osd_351, \
                         osd_353 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_586[k] = pa_x[k] * nsf0_586[k]
                   - f_4 * pc_x[k] * nsf1_586[k];

        t_587[k] = f_12 * nsd_285[k]
                   + f_3 * pc_z[k] * osd_351[k];

        t_588[k] = f_13 * nsd_293[k]
                   + f_3 * pc_y[k] * osd_353[k];

        t_589[k] = pa_x[k] * nsf0_589[k]
                   - f_4 * pc_x[k] * nsf1_589[k];
    }

#pragma omp simd aligned(t_590, t_591, t_592, t_593, pa_x, pc_x, pc_y, pc_z, nsf0_590, \
                         nsd_288, nsd_294, nsd_354, nsd_357, nsf1_590, osd_354, \
                         osd_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_590[k] = pa_x[k] * nsf0_590[k]
                   + f_12 * nsd_354[k]
                   - f_4 * pc_x[k] * nsf1_590[k];

        t_591[k] = f_15 * nsd_294[k]
                   + f_3 * pc_y[k] * osd_354[k];

        t_592[k] = f_14 * nsd_288[k]
                   + f_3 * pc_z[k] * osd_354[k];

        t_593[k] = f_5 * nsd_357[k]
                   + f_3 * pc_x[k] * osd_357[k];
    }

#pragma omp simd aligned(t_594, t_595, t_596, t_597, pa_x, pc_x, pc_z, nsf0_596, nsd_291, \
                         nsd_358, nsd_359, nsf1_596, osd_357, osd_358, \
                         osd_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_594[k] = f_5 * nsd_358[k]
                   + f_3 * pc_x[k] * osd_358[k];

        t_595[k] = f_5 * nsd_359[k]
                   + f_3 * pc_x[k] * osd_359[k];

        t_596[k] = pa_x[k] * nsf0_596[k]
                   - f_4 * pc_x[k] * nsf1_596[k];

        t_597[k] = f_14 * nsd_291[k]
                   + f_3 * pc_z[k] * osd_357[k];
    }

#pragma omp simd aligned(t_598, t_599, t_600, t_601, pa_x, pc_x, pc_y, nsf0_599, nsf0_600, \
                         nsd_299, nsd_300, nsd_360, nsf1_599, nsf1_600, osd_359, \
                         osd_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_598[k] = f_15 * nsd_299[k]
                   + f_3 * pc_y[k] * osd_359[k];

        t_599[k] = pa_x[k] * nsf0_599[k]
                   - f_4 * pc_x[k] * nsf1_599[k];

        t_600[k] = pa_x[k] * nsf0_600[k]
                   + f_12 * nsd_360[k]
                   - f_4 * pc_x[k] * nsf1_600[k];

        t_601[k] = f_16 * nsd_300[k]
                   + f_3 * pc_y[k] * osd_360[k];
    }

#pragma omp simd aligned(t_602, t_603, t_604, t_605, pc_x, pc_z, nsd_294, nsd_363, nsd_364, \
                         nsd_365, osd_360, osd_363, osd_364, osd_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_602[k] = f_16 * nsd_294[k]
                   + f_3 * pc_z[k] * osd_360[k];

        t_603[k] = f_5 * nsd_363[k]
                   + f_3 * pc_x[k] * osd_363[k];

        t_604[k] = f_5 * nsd_364[k]
                   + f_3 * pc_x[k] * osd_364[k];

        t_605[k] = f_5 * nsd_365[k]
                   + f_3 * pc_x[k] * osd_365[k];
    }

#pragma omp simd aligned(t_606, t_607, t_608, t_609, pa_x, pc_x, pc_y, pc_z, nsf0_606, \
                         nsf0_609, nsd_297, nsd_305, nsf1_606, nsf1_609, osd_363, \
                         osd_365 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_606[k] = pa_x[k] * nsf0_606[k]
                   - f_4 * pc_x[k] * nsf1_606[k];

        t_607[k] = f_16 * nsd_297[k]
                   + f_3 * pc_z[k] * osd_363[k];

        t_608[k] = f_16 * nsd_305[k]
                   + f_3 * pc_y[k] * osd_365[k];

        t_609[k] = pa_x[k] * nsf0_609[k]
                   - f_4 * pc_x[k] * nsf1_609[k];
    }

#pragma omp simd aligned(t_610, t_611, t_612, t_613, pa_x, pc_x, pc_y, pc_z, nsf0_610, \
                         nsd_300, nsd_306, nsd_366, nsd_369, nsf1_610, osd_366, \
                         osd_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_610[k] = pa_x[k] * nsf0_610[k]
                   + f_12 * nsd_366[k]
                   - f_4 * pc_x[k] * nsf1_610[k];

        t_611[k] = f_14 * nsd_306[k]
                   + f_3 * pc_y[k] * osd_366[k];

        t_612[k] = f_15 * nsd_300[k]
                   + f_3 * pc_z[k] * osd_366[k];

        t_613[k] = f_5 * nsd_369[k]
                   + f_3 * pc_x[k] * osd_369[k];
    }

#pragma omp simd aligned(t_614, t_615, t_616, t_617, pa_x, pc_x, pc_z, nsf0_616, nsd_303, \
                         nsd_370, nsd_371, nsf1_616, osd_369, osd_370, \
                         osd_371 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_614[k] = f_5 * nsd_370[k]
                   + f_3 * pc_x[k] * osd_370[k];

        t_615[k] = f_5 * nsd_371[k]
                   + f_3 * pc_x[k] * osd_371[k];

        t_616[k] = pa_x[k] * nsf0_616[k]
                   - f_4 * pc_x[k] * nsf1_616[k];

        t_617[k] = f_15 * nsd_303[k]
                   + f_3 * pc_z[k] * osd_369[k];
    }

#pragma omp simd aligned(t_618, t_619, t_620, t_621, pa_x, pc_x, pc_y, nsf0_619, nsf0_620, \
                         nsd_311, nsd_312, nsd_372, nsf1_619, nsf1_620, osd_371, \
                         osd_372 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_618[k] = f_14 * nsd_311[k]
                   + f_3 * pc_y[k] * osd_371[k];

        t_619[k] = pa_x[k] * nsf0_619[k]
                   - f_4 * pc_x[k] * nsf1_619[k];

        t_620[k] = pa_x[k] * nsf0_620[k]
                   + f_12 * nsd_372[k]
                   - f_4 * pc_x[k] * nsf1_620[k];

        t_621[k] = f_12 * nsd_312[k]
                   + f_3 * pc_y[k] * osd_372[k];
    }

#pragma omp simd aligned(t_622, t_623, t_624, t_625, pc_x, pc_z, nsd_306, nsd_375, nsd_376, \
                         nsd_377, osd_372, osd_375, osd_376, osd_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_622[k] = f_13 * nsd_306[k]
                   + f_3 * pc_z[k] * osd_372[k];

        t_623[k] = f_5 * nsd_375[k]
                   + f_3 * pc_x[k] * osd_375[k];

        t_624[k] = f_5 * nsd_376[k]
                   + f_3 * pc_x[k] * osd_376[k];

        t_625[k] = f_5 * nsd_377[k]
                   + f_3 * pc_x[k] * osd_377[k];
    }
}

static auto
compute_prim_osf_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsf0,
                                                          const size_t nsd, const size_t nsf1,
                                                          const size_t osp0, const size_t osp1,
                                                          const size_t osd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 5.0 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_9 = 4.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_11 = 4.0 / q;
    const auto f_12 = 1.5 / q;
    const auto f_13 = 3.5 / q;
    const auto f_14 = 2.0 / q;
    const auto f_15 = 3.0 / q;
    const auto f_16 = 2.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsf0_540 = buffer.data(nsf0 + 540);
    const auto *nsf0_550 = buffer.data(nsf0 + 550);
    const auto *nsf0_551 = buffer.data(nsf0 + 551);
    const auto *nsf0_556 = buffer.data(nsf0 + 556);
    const auto *nsf0_626 = buffer.data(nsf0 + 626);
    const auto *nsf0_629 = buffer.data(nsf0 + 629);
    const auto *nsf0_630 = buffer.data(nsf0 + 630);
    const auto *nsf0_636 = buffer.data(nsf0 + 636);
    const auto *nsf0_639 = buffer.data(nsf0 + 639);
    const auto *nsf0_646 = buffer.data(nsf0 + 646);
    const auto *nsf0_649 = buffer.data(nsf0 + 649);
    const auto *nsf0_650 = buffer.data(nsf0 + 650);
    const auto *nsf0_656 = buffer.data(nsf0 + 656);
    const auto *nsf0_657 = buffer.data(nsf0 + 657);
    const auto *nsf0_659 = buffer.data(nsf0 + 659);

    const auto *nsd_309 = buffer.data(nsd + 309);
    const auto *nsd_312 = buffer.data(nsd + 312);
    const auto *nsd_315 = buffer.data(nsd + 315);
    const auto *nsd_317 = buffer.data(nsd + 317);
    const auto *nsd_318 = buffer.data(nsd + 318);
    const auto *nsd_321 = buffer.data(nsd + 321);
    const auto *nsd_323 = buffer.data(nsd + 323);
    const auto *nsd_324 = buffer.data(nsd + 324);
    const auto *nsd_329 = buffer.data(nsd + 329);
    const auto *nsd_333 = buffer.data(nsd + 333);
    const auto *nsd_335 = buffer.data(nsd + 335);
    const auto *nsd_339 = buffer.data(nsd + 339);
    const auto *nsd_341 = buffer.data(nsd + 341);
    const auto *nsd_345 = buffer.data(nsd + 345);
    const auto *nsd_347 = buffer.data(nsd + 347);
    const auto *nsd_351 = buffer.data(nsd + 351);
    const auto *nsd_353 = buffer.data(nsd + 353);
    const auto *nsd_357 = buffer.data(nsd + 357);
    const auto *nsd_359 = buffer.data(nsd + 359);
    const auto *nsd_363 = buffer.data(nsd + 363);
    const auto *nsd_365 = buffer.data(nsd + 365);
    const auto *nsd_369 = buffer.data(nsd + 369);
    const auto *nsd_371 = buffer.data(nsd + 371);
    const auto *nsd_375 = buffer.data(nsd + 375);
    const auto *nsd_377 = buffer.data(nsd + 377);
    const auto *nsd_378 = buffer.data(nsd + 378);
    const auto *nsd_381 = buffer.data(nsd + 381);
    const auto *nsd_382 = buffer.data(nsd + 382);
    const auto *nsd_383 = buffer.data(nsd + 383);
    const auto *nsd_387 = buffer.data(nsd + 387);
    const auto *nsd_388 = buffer.data(nsd + 388);
    const auto *nsd_389 = buffer.data(nsd + 389);
    const auto *nsd_390 = buffer.data(nsd + 390);
    const auto *nsd_393 = buffer.data(nsd + 393);
    const auto *nsd_395 = buffer.data(nsd + 395);

    const auto *nsf1_540 = buffer.data(nsf1 + 540);
    const auto *nsf1_550 = buffer.data(nsf1 + 550);
    const auto *nsf1_551 = buffer.data(nsf1 + 551);
    const auto *nsf1_556 = buffer.data(nsf1 + 556);
    const auto *nsf1_626 = buffer.data(nsf1 + 626);
    const auto *nsf1_629 = buffer.data(nsf1 + 629);
    const auto *nsf1_630 = buffer.data(nsf1 + 630);
    const auto *nsf1_636 = buffer.data(nsf1 + 636);
    const auto *nsf1_639 = buffer.data(nsf1 + 639);
    const auto *nsf1_646 = buffer.data(nsf1 + 646);
    const auto *nsf1_649 = buffer.data(nsf1 + 649);
    const auto *nsf1_650 = buffer.data(nsf1 + 650);
    const auto *nsf1_656 = buffer.data(nsf1 + 656);
    const auto *nsf1_657 = buffer.data(nsf1 + 657);
    const auto *nsf1_659 = buffer.data(nsf1 + 659);

    const auto *osp0_198 = buffer.data(osp0 + 198);
    const auto *osp0_199 = buffer.data(osp0 + 199);
    const auto *osp0_200 = buffer.data(osp0 + 200);
    const auto *osp0_203 = buffer.data(osp0 + 203);
    const auto *osp0_204 = buffer.data(osp0 + 204);
    const auto *osp0_205 = buffer.data(osp0 + 205);
    const auto *osp0_206 = buffer.data(osp0 + 206);
    const auto *osp0_207 = buffer.data(osp0 + 207);
    const auto *osp0_208 = buffer.data(osp0 + 208);
    const auto *osp0_209 = buffer.data(osp0 + 209);
    const auto *osp0_210 = buffer.data(osp0 + 210);
    const auto *osp0_211 = buffer.data(osp0 + 211);
    const auto *osp0_212 = buffer.data(osp0 + 212);
    const auto *osp0_213 = buffer.data(osp0 + 213);
    const auto *osp0_214 = buffer.data(osp0 + 214);
    const auto *osp0_215 = buffer.data(osp0 + 215);
    const auto *osp0_216 = buffer.data(osp0 + 216);
    const auto *osp0_217 = buffer.data(osp0 + 217);
    const auto *osp0_218 = buffer.data(osp0 + 218);
    const auto *osp0_219 = buffer.data(osp0 + 219);
    const auto *osp0_220 = buffer.data(osp0 + 220);
    const auto *osp0_221 = buffer.data(osp0 + 221);
    const auto *osp0_222 = buffer.data(osp0 + 222);
    const auto *osp0_223 = buffer.data(osp0 + 223);
    const auto *osp0_224 = buffer.data(osp0 + 224);
    const auto *osp0_225 = buffer.data(osp0 + 225);
    const auto *osp0_226 = buffer.data(osp0 + 226);
    const auto *osp0_227 = buffer.data(osp0 + 227);

    const auto *osp1_198 = buffer.data(osp1 + 198);
    const auto *osp1_199 = buffer.data(osp1 + 199);
    const auto *osp1_200 = buffer.data(osp1 + 200);
    const auto *osp1_203 = buffer.data(osp1 + 203);
    const auto *osp1_204 = buffer.data(osp1 + 204);
    const auto *osp1_205 = buffer.data(osp1 + 205);
    const auto *osp1_206 = buffer.data(osp1 + 206);
    const auto *osp1_207 = buffer.data(osp1 + 207);
    const auto *osp1_208 = buffer.data(osp1 + 208);
    const auto *osp1_209 = buffer.data(osp1 + 209);
    const auto *osp1_210 = buffer.data(osp1 + 210);
    const auto *osp1_211 = buffer.data(osp1 + 211);
    const auto *osp1_212 = buffer.data(osp1 + 212);
    const auto *osp1_213 = buffer.data(osp1 + 213);
    const auto *osp1_214 = buffer.data(osp1 + 214);
    const auto *osp1_215 = buffer.data(osp1 + 215);
    const auto *osp1_216 = buffer.data(osp1 + 216);
    const auto *osp1_217 = buffer.data(osp1 + 217);
    const auto *osp1_218 = buffer.data(osp1 + 218);
    const auto *osp1_219 = buffer.data(osp1 + 219);
    const auto *osp1_220 = buffer.data(osp1 + 220);
    const auto *osp1_221 = buffer.data(osp1 + 221);
    const auto *osp1_222 = buffer.data(osp1 + 222);
    const auto *osp1_223 = buffer.data(osp1 + 223);
    const auto *osp1_224 = buffer.data(osp1 + 224);
    const auto *osp1_225 = buffer.data(osp1 + 225);
    const auto *osp1_226 = buffer.data(osp1 + 226);
    const auto *osp1_227 = buffer.data(osp1 + 227);

    const auto *osd_375 = buffer.data(osd + 375);
    const auto *osd_377 = buffer.data(osd + 377);
    const auto *osd_378 = buffer.data(osd + 378);
    const auto *osd_381 = buffer.data(osd + 381);
    const auto *osd_382 = buffer.data(osd + 382);
    const auto *osd_383 = buffer.data(osd + 383);
    const auto *osd_384 = buffer.data(osd + 384);
    const auto *osd_387 = buffer.data(osd + 387);
    const auto *osd_388 = buffer.data(osd + 388);
    const auto *osd_389 = buffer.data(osd + 389);
    const auto *osd_390 = buffer.data(osd + 390);
    const auto *osd_392 = buffer.data(osd + 392);
    const auto *osd_393 = buffer.data(osd + 393);
    const auto *osd_395 = buffer.data(osd + 395);
    const auto *osd_396 = buffer.data(osd + 396);
    const auto *osd_397 = buffer.data(osd + 397);
    const auto *osd_399 = buffer.data(osd + 399);
    const auto *osd_400 = buffer.data(osd + 400);
    const auto *osd_401 = buffer.data(osd + 401);
    const auto *osd_404 = buffer.data(osd + 404);
    const auto *osd_405 = buffer.data(osd + 405);
    const auto *osd_406 = buffer.data(osd + 406);
    const auto *osd_407 = buffer.data(osd + 407);
    const auto *osd_408 = buffer.data(osd + 408);
    const auto *osd_409 = buffer.data(osd + 409);
    const auto *osd_410 = buffer.data(osd + 410);
    const auto *osd_411 = buffer.data(osd + 411);
    const auto *osd_412 = buffer.data(osd + 412);
    const auto *osd_413 = buffer.data(osd + 413);
    const auto *osd_414 = buffer.data(osd + 414);
    const auto *osd_415 = buffer.data(osd + 415);
    const auto *osd_416 = buffer.data(osd + 416);
    const auto *osd_417 = buffer.data(osd + 417);
    const auto *osd_418 = buffer.data(osd + 418);
    const auto *osd_419 = buffer.data(osd + 419);
    const auto *osd_420 = buffer.data(osd + 420);
    const auto *osd_421 = buffer.data(osd + 421);
    const auto *osd_422 = buffer.data(osd + 422);
    const auto *osd_423 = buffer.data(osd + 423);
    const auto *osd_424 = buffer.data(osd + 424);
    const auto *osd_425 = buffer.data(osd + 425);
    const auto *osd_426 = buffer.data(osd + 426);
    const auto *osd_427 = buffer.data(osd + 427);
    const auto *osd_428 = buffer.data(osd + 428);
    const auto *osd_429 = buffer.data(osd + 429);
    const auto *osd_430 = buffer.data(osd + 430);
    const auto *osd_431 = buffer.data(osd + 431);
    const auto *osd_432 = buffer.data(osd + 432);
    const auto *osd_433 = buffer.data(osd + 433);
    const auto *osd_434 = buffer.data(osd + 434);
    const auto *osd_435 = buffer.data(osd + 435);
    const auto *osd_436 = buffer.data(osd + 436);
    const auto *osd_437 = buffer.data(osd + 437);
    const auto *osd_438 = buffer.data(osd + 438);
    const auto *osd_439 = buffer.data(osd + 439);
    const auto *osd_440 = buffer.data(osd + 440);
    const auto *osd_441 = buffer.data(osd + 441);
    const auto *osd_442 = buffer.data(osd + 442);
    const auto *osd_443 = buffer.data(osd + 443);
    const auto *osd_444 = buffer.data(osd + 444);
    const auto *osd_445 = buffer.data(osd + 445);
    const auto *osd_446 = buffer.data(osd + 446);
    const auto *osd_447 = buffer.data(osd + 447);
    const auto *osd_448 = buffer.data(osd + 448);
    const auto *osd_449 = buffer.data(osd + 449);
    const auto *osd_450 = buffer.data(osd + 450);
    const auto *osd_451 = buffer.data(osd + 451);
    const auto *osd_452 = buffer.data(osd + 452);
    const auto *osd_453 = buffer.data(osd + 453);
    const auto *osd_454 = buffer.data(osd + 454);
    const auto *osd_455 = buffer.data(osd + 455);

#pragma omp simd aligned(t_626, t_627, t_628, t_629, pa_x, pc_x, pc_y, pc_z, nsf0_626, \
                         nsf0_629, nsd_309, nsd_317, nsf1_626, nsf1_629, osd_375, \
                         osd_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = pa_x[k] * nsf0_626[k]
                   - f_4 * pc_x[k] * nsf1_626[k];

        t_627[k] = f_13 * nsd_309[k]
                   + f_3 * pc_z[k] * osd_375[k];

        t_628[k] = f_12 * nsd_317[k]
                   + f_3 * pc_y[k] * osd_377[k];

        t_629[k] = pa_x[k] * nsf0_629[k]
                   - f_4 * pc_x[k] * nsf1_629[k];
    }

#pragma omp simd aligned(t_630, t_631, t_632, t_633, pa_x, pc_x, pc_y, pc_z, nsf0_630, \
                         nsd_312, nsd_318, nsd_378, nsd_381, nsf1_630, osd_378, \
                         osd_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_630[k] = pa_x[k] * nsf0_630[k]
                   + f_12 * nsd_378[k]
                   - f_4 * pc_x[k] * nsf1_630[k];

        t_631[k] = f_10 * nsd_318[k]
                   + f_3 * pc_y[k] * osd_378[k];

        t_632[k] = f_11 * nsd_312[k]
                   + f_3 * pc_z[k] * osd_378[k];

        t_633[k] = f_5 * nsd_381[k]
                   + f_3 * pc_x[k] * osd_381[k];
    }

#pragma omp simd aligned(t_634, t_635, t_636, t_637, pa_x, pc_x, pc_z, nsf0_636, nsd_315, \
                         nsd_382, nsd_383, nsf1_636, osd_381, osd_382, \
                         osd_383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_634[k] = f_5 * nsd_382[k]
                   + f_3 * pc_x[k] * osd_382[k];

        t_635[k] = f_5 * nsd_383[k]
                   + f_3 * pc_x[k] * osd_383[k];

        t_636[k] = pa_x[k] * nsf0_636[k]
                   - f_4 * pc_x[k] * nsf1_636[k];

        t_637[k] = f_11 * nsd_315[k]
                   + f_3 * pc_z[k] * osd_381[k];
    }

#pragma omp simd aligned(t_638, t_639, t_640, t_641, pa_x, pa_y, pc_x, pc_y, nsf0_540, \
                         nsf0_639, nsd_323, nsd_324, nsf1_540, nsf1_639, osd_383, \
                         osd_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_638[k] = f_10 * nsd_323[k]
                   + f_3 * pc_y[k] * osd_383[k];

        t_639[k] = pa_x[k] * nsf0_639[k]
                   - f_4 * pc_x[k] * nsf1_639[k];

        t_640[k] = pa_y[k] * nsf0_540[k]
                   - f_4 * pc_y[k] * nsf1_540[k];

        t_641[k] = f_5 * nsd_324[k]
                   + f_3 * pc_y[k] * osd_384[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, t_645, pc_x, pc_z, nsd_318, nsd_387, nsd_388, \
                         nsd_389, osd_384, osd_387, osd_388, osd_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = f_9 * nsd_318[k]
                   + f_3 * pc_z[k] * osd_384[k];

        t_643[k] = f_5 * nsd_387[k]
                   + f_3 * pc_x[k] * osd_387[k];

        t_644[k] = f_5 * nsd_388[k]
                   + f_3 * pc_x[k] * osd_388[k];

        t_645[k] = f_5 * nsd_389[k]
                   + f_3 * pc_x[k] * osd_389[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, t_649, pa_x, pc_x, pc_y, pc_z, nsf0_646, \
                         nsf0_649, nsd_321, nsd_329, nsf1_646, nsf1_649, osd_387, \
                         osd_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = pa_x[k] * nsf0_646[k]
                   - f_4 * pc_x[k] * nsf1_646[k];

        t_647[k] = f_9 * nsd_321[k]
                   + f_3 * pc_z[k] * osd_387[k];

        t_648[k] = f_5 * nsd_329[k]
                   + f_3 * pc_y[k] * osd_389[k];

        t_649[k] = pa_x[k] * nsf0_649[k]
                   - f_4 * pc_x[k] * nsf1_649[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, t_653, pa_x, pc_x, pc_y, pc_z, nsf0_650, \
                         nsd_324, nsd_390, nsd_393, nsf1_650, osd_390, \
                         osd_393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = pa_x[k] * nsf0_650[k]
                   + f_12 * nsd_390[k]
                   - f_4 * pc_x[k] * nsf1_650[k];

        t_651[k] = f_3 * pc_y[k] * osd_390[k];

        t_652[k] = f_6 * nsd_324[k]
                   + f_3 * pc_z[k] * osd_390[k];

        t_653[k] = f_5 * nsd_393[k]
                   + f_3 * pc_x[k] * osd_393[k];
    }

#pragma omp simd aligned(t_654, t_655, t_656, t_657, t_658, pa_x, pc_x, pc_y, nsf0_656, \
                         nsf0_657, nsd_395, nsf1_656, nsf1_657, osd_392, \
                         osd_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_654[k] = f_3 * pc_y[k] * osd_392[k];

        t_655[k] = f_5 * nsd_395[k]
                   + f_3 * pc_x[k] * osd_395[k];

        t_656[k] = pa_x[k] * nsf0_656[k]
                   - f_4 * pc_x[k] * nsf1_656[k];

        t_657[k] = pa_x[k] * nsf0_657[k]
                   - f_4 * pc_x[k] * nsf1_657[k];

        t_658[k] = f_3 * pc_y[k] * osd_395[k];
    }

#pragma omp simd aligned(t_659, t_660, t_661, t_662, pa_x, pc_x, pc_z, nsf0_659, nsf1_659, \
                         osp0_198, osp0_199, osp1_198, osp1_199, osd_396, \
                         osd_397 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_659[k] = pa_x[k] * nsf0_659[k]
                   - f_4 * pc_x[k] * nsf1_659[k];

        t_660[k] = f_1 * osp0_198[k]
                   - f_2 * osp1_198[k]
                   + f_3 * pc_x[k] * osd_396[k];

        t_661[k] = f_7 * osp0_199[k]
                   - f_8 * osp1_199[k]
                   + f_3 * pc_x[k] * osd_397[k];

        t_662[k] = f_3 * pc_z[k] * osd_396[k];
    }

#pragma omp simd aligned(t_663, t_664, t_665, t_666, t_667, t_668, pc_x, pc_y, pc_z, nsd_333, \
                         nsd_335, osp0_199, osp1_199, osd_399, osd_400, \
                         osd_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_663[k] = f_3 * pc_x[k] * osd_399[k];

        t_664[k] = f_3 * pc_x[k] * osd_400[k];

        t_665[k] = f_3 * pc_x[k] * osd_401[k];

        t_666[k] = f_0 * nsd_333[k]
                   + f_1 * osp0_199[k]
                   - f_2 * osp1_199[k]
                   + f_3 * pc_y[k] * osd_399[k];

        t_667[k] = f_3 * pc_z[k] * osd_399[k];

        t_668[k] = f_0 * nsd_335[k]
                   + f_3 * pc_y[k] * osd_401[k];
    }

#pragma omp simd aligned(t_669, t_670, t_671, pa_z, pc_z, nsf0_550, nsf0_551, nsf1_550, \
                         nsf1_551, osp0_200, osp1_200, osd_401 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_669[k] = f_1 * osp0_200[k]
                   - f_2 * osp1_200[k]
                   + f_3 * pc_z[k] * osd_401[k];

        t_670[k] = pa_z[k] * nsf0_550[k]
                   - f_4 * pc_z[k] * nsf1_550[k];

        t_671[k] = pa_z[k] * nsf0_551[k]
                   - f_4 * pc_z[k] * nsf1_551[k];
    }

#pragma omp simd aligned(t_672, t_673, t_674, t_675, t_676, pa_z, pc_x, pc_z, nsf0_556, \
                         nsf1_556, osp0_203, osp1_203, osd_404, osd_405, osd_406, \
                         osd_407 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_672[k] = f_7 * osp0_203[k]
                   - f_8 * osp1_203[k]
                   + f_3 * pc_x[k] * osd_404[k];

        t_673[k] = f_3 * pc_x[k] * osd_405[k];

        t_674[k] = f_3 * pc_x[k] * osd_406[k];

        t_675[k] = f_3 * pc_x[k] * osd_407[k];

        t_676[k] = pa_z[k] * nsf0_556[k]
                   - f_4 * pc_z[k] * nsf1_556[k];
    }

#pragma omp simd aligned(t_677, t_678, t_679, pc_y, pc_z, nsd_333, nsd_335, nsd_341, osp0_203, \
                         osp1_203, osd_405, osd_407 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_677[k] = f_5 * nsd_333[k]
                   + f_3 * pc_z[k] * osd_405[k];

        t_678[k] = f_6 * nsd_341[k]
                   + f_3 * pc_y[k] * osd_407[k];

        t_679[k] = f_5 * nsd_335[k]
                   + f_1 * osp0_203[k]
                   - f_2 * osp1_203[k]
                   + f_3 * pc_z[k] * osd_407[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, t_683, pc_x, osp0_204, osp0_205, osp0_206, \
                         osp1_204, osp1_205, osp1_206, osd_408, osd_409, osd_410, \
                         osd_411 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_1 * osp0_204[k]
                   - f_2 * osp1_204[k]
                   + f_3 * pc_x[k] * osd_408[k];

        t_681[k] = f_7 * osp0_205[k]
                   - f_8 * osp1_205[k]
                   + f_3 * pc_x[k] * osd_409[k];

        t_682[k] = f_7 * osp0_206[k]
                   - f_8 * osp1_206[k]
                   + f_3 * pc_x[k] * osd_410[k];

        t_683[k] = f_3 * pc_x[k] * osd_411[k];
    }

#pragma omp simd aligned(t_684, t_685, t_686, t_687, t_688, pc_x, pc_y, pc_z, nsd_339, \
                         nsd_345, nsd_347, osp0_205, osp1_205, osd_411, osd_412, \
                         osd_413 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_684[k] = f_3 * pc_x[k] * osd_412[k];

        t_685[k] = f_3 * pc_x[k] * osd_413[k];

        t_686[k] = f_9 * nsd_345[k]
                   + f_1 * osp0_205[k]
                   - f_2 * osp1_205[k]
                   + f_3 * pc_y[k] * osd_411[k];

        t_687[k] = f_10 * nsd_339[k]
                   + f_3 * pc_z[k] * osd_411[k];

        t_688[k] = f_9 * nsd_347[k]
                   + f_3 * pc_y[k] * osd_413[k];
    }

#pragma omp simd aligned(t_689, t_690, t_691, pc_x, pc_z, nsd_341, osp0_206, osp0_207, \
                         osp0_208, osp1_206, osp1_207, osp1_208, osd_413, osd_414, \
                         osd_415 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_689[k] = f_10 * nsd_341[k]
                   + f_1 * osp0_206[k]
                   - f_2 * osp1_206[k]
                   + f_3 * pc_z[k] * osd_413[k];

        t_690[k] = f_1 * osp0_207[k]
                   - f_2 * osp1_207[k]
                   + f_3 * pc_x[k] * osd_414[k];

        t_691[k] = f_7 * osp0_208[k]
                   - f_8 * osp1_208[k]
                   + f_3 * pc_x[k] * osd_415[k];
    }

#pragma omp simd aligned(t_692, t_693, t_694, t_695, t_696, pc_x, pc_y, nsd_351, osp0_208, \
                         osp0_209, osp1_208, osp1_209, osd_416, osd_417, osd_418, \
                         osd_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_692[k] = f_7 * osp0_209[k]
                   - f_8 * osp1_209[k]
                   + f_3 * pc_x[k] * osd_416[k];

        t_693[k] = f_3 * pc_x[k] * osd_417[k];

        t_694[k] = f_3 * pc_x[k] * osd_418[k];

        t_695[k] = f_3 * pc_x[k] * osd_419[k];

        t_696[k] = f_11 * nsd_351[k]
                   + f_1 * osp0_208[k]
                   - f_2 * osp1_208[k]
                   + f_3 * pc_y[k] * osd_417[k];
    }

#pragma omp simd aligned(t_697, t_698, t_699, pc_y, pc_z, nsd_345, nsd_347, nsd_353, osp0_209, \
                         osp1_209, osd_417, osd_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_697[k] = f_12 * nsd_345[k]
                   + f_3 * pc_z[k] * osd_417[k];

        t_698[k] = f_11 * nsd_353[k]
                   + f_3 * pc_y[k] * osd_419[k];

        t_699[k] = f_12 * nsd_347[k]
                   + f_1 * osp0_209[k]
                   - f_2 * osp1_209[k]
                   + f_3 * pc_z[k] * osd_419[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, t_703, pc_x, osp0_210, osp0_211, osp0_212, \
                         osp1_210, osp1_211, osp1_212, osd_420, osd_421, osd_422, \
                         osd_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = f_1 * osp0_210[k]
                   - f_2 * osp1_210[k]
                   + f_3 * pc_x[k] * osd_420[k];

        t_701[k] = f_7 * osp0_211[k]
                   - f_8 * osp1_211[k]
                   + f_3 * pc_x[k] * osd_421[k];

        t_702[k] = f_7 * osp0_212[k]
                   - f_8 * osp1_212[k]
                   + f_3 * pc_x[k] * osd_422[k];

        t_703[k] = f_3 * pc_x[k] * osd_423[k];
    }

#pragma omp simd aligned(t_704, t_705, t_706, t_707, t_708, pc_x, pc_y, pc_z, nsd_351, \
                         nsd_357, nsd_359, osp0_211, osp1_211, osd_423, osd_424, \
                         osd_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_704[k] = f_3 * pc_x[k] * osd_424[k];

        t_705[k] = f_3 * pc_x[k] * osd_425[k];

        t_706[k] = f_13 * nsd_357[k]
                   + f_1 * osp0_211[k]
                   - f_2 * osp1_211[k]
                   + f_3 * pc_y[k] * osd_423[k];

        t_707[k] = f_14 * nsd_351[k]
                   + f_3 * pc_z[k] * osd_423[k];

        t_708[k] = f_13 * nsd_359[k]
                   + f_3 * pc_y[k] * osd_425[k];
    }

#pragma omp simd aligned(t_709, t_710, t_711, pc_x, pc_z, nsd_353, osp0_212, osp0_213, \
                         osp0_214, osp1_212, osp1_213, osp1_214, osd_425, osd_426, \
                         osd_427 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_709[k] = f_14 * nsd_353[k]
                   + f_1 * osp0_212[k]
                   - f_2 * osp1_212[k]
                   + f_3 * pc_z[k] * osd_425[k];

        t_710[k] = f_1 * osp0_213[k]
                   - f_2 * osp1_213[k]
                   + f_3 * pc_x[k] * osd_426[k];

        t_711[k] = f_7 * osp0_214[k]
                   - f_8 * osp1_214[k]
                   + f_3 * pc_x[k] * osd_427[k];
    }

#pragma omp simd aligned(t_712, t_713, t_714, t_715, t_716, pc_x, pc_y, nsd_363, osp0_214, \
                         osp0_215, osp1_214, osp1_215, osd_428, osd_429, osd_430, \
                         osd_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = f_7 * osp0_215[k]
                   - f_8 * osp1_215[k]
                   + f_3 * pc_x[k] * osd_428[k];

        t_713[k] = f_3 * pc_x[k] * osd_429[k];

        t_714[k] = f_3 * pc_x[k] * osd_430[k];

        t_715[k] = f_3 * pc_x[k] * osd_431[k];

        t_716[k] = f_15 * nsd_363[k]
                   + f_1 * osp0_214[k]
                   - f_2 * osp1_214[k]
                   + f_3 * pc_y[k] * osd_429[k];
    }

#pragma omp simd aligned(t_717, t_718, t_719, pc_y, pc_z, nsd_357, nsd_359, nsd_365, osp0_215, \
                         osp1_215, osd_429, osd_431 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_717[k] = f_16 * nsd_357[k]
                   + f_3 * pc_z[k] * osd_429[k];

        t_718[k] = f_15 * nsd_365[k]
                   + f_3 * pc_y[k] * osd_431[k];

        t_719[k] = f_16 * nsd_359[k]
                   + f_1 * osp0_215[k]
                   - f_2 * osp1_215[k]
                   + f_3 * pc_z[k] * osd_431[k];
    }

#pragma omp simd aligned(t_720, t_721, t_722, t_723, pc_x, osp0_216, osp0_217, osp0_218, \
                         osp1_216, osp1_217, osp1_218, osd_432, osd_433, osd_434, \
                         osd_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = f_1 * osp0_216[k]
                   - f_2 * osp1_216[k]
                   + f_3 * pc_x[k] * osd_432[k];

        t_721[k] = f_7 * osp0_217[k]
                   - f_8 * osp1_217[k]
                   + f_3 * pc_x[k] * osd_433[k];

        t_722[k] = f_7 * osp0_218[k]
                   - f_8 * osp1_218[k]
                   + f_3 * pc_x[k] * osd_434[k];

        t_723[k] = f_3 * pc_x[k] * osd_435[k];
    }

#pragma omp simd aligned(t_724, t_725, t_726, t_727, t_728, pc_x, pc_y, pc_z, nsd_363, \
                         nsd_369, nsd_371, osp0_217, osp1_217, osd_435, osd_436, \
                         osd_437 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_724[k] = f_3 * pc_x[k] * osd_436[k];

        t_725[k] = f_3 * pc_x[k] * osd_437[k];

        t_726[k] = f_16 * nsd_369[k]
                   + f_1 * osp0_217[k]
                   - f_2 * osp1_217[k]
                   + f_3 * pc_y[k] * osd_435[k];

        t_727[k] = f_15 * nsd_363[k]
                   + f_3 * pc_z[k] * osd_435[k];

        t_728[k] = f_16 * nsd_371[k]
                   + f_3 * pc_y[k] * osd_437[k];
    }

#pragma omp simd aligned(t_729, t_730, t_731, pc_x, pc_z, nsd_365, osp0_218, osp0_219, \
                         osp0_220, osp1_218, osp1_219, osp1_220, osd_437, osd_438, \
                         osd_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_729[k] = f_15 * nsd_365[k]
                   + f_1 * osp0_218[k]
                   - f_2 * osp1_218[k]
                   + f_3 * pc_z[k] * osd_437[k];

        t_730[k] = f_1 * osp0_219[k]
                   - f_2 * osp1_219[k]
                   + f_3 * pc_x[k] * osd_438[k];

        t_731[k] = f_7 * osp0_220[k]
                   - f_8 * osp1_220[k]
                   + f_3 * pc_x[k] * osd_439[k];
    }

#pragma omp simd aligned(t_732, t_733, t_734, t_735, t_736, pc_x, pc_y, nsd_375, osp0_220, \
                         osp0_221, osp1_220, osp1_221, osd_440, osd_441, osd_442, \
                         osd_443 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_732[k] = f_7 * osp0_221[k]
                   - f_8 * osp1_221[k]
                   + f_3 * pc_x[k] * osd_440[k];

        t_733[k] = f_3 * pc_x[k] * osd_441[k];

        t_734[k] = f_3 * pc_x[k] * osd_442[k];

        t_735[k] = f_3 * pc_x[k] * osd_443[k];

        t_736[k] = f_14 * nsd_375[k]
                   + f_1 * osp0_220[k]
                   - f_2 * osp1_220[k]
                   + f_3 * pc_y[k] * osd_441[k];
    }

#pragma omp simd aligned(t_737, t_738, t_739, pc_y, pc_z, nsd_369, nsd_371, nsd_377, osp0_221, \
                         osp1_221, osd_441, osd_443 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_737[k] = f_13 * nsd_369[k]
                   + f_3 * pc_z[k] * osd_441[k];

        t_738[k] = f_14 * nsd_377[k]
                   + f_3 * pc_y[k] * osd_443[k];

        t_739[k] = f_13 * nsd_371[k]
                   + f_1 * osp0_221[k]
                   - f_2 * osp1_221[k]
                   + f_3 * pc_z[k] * osd_443[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, t_743, pc_x, osp0_222, osp0_223, osp0_224, \
                         osp1_222, osp1_223, osp1_224, osd_444, osd_445, osd_446, \
                         osd_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = f_1 * osp0_222[k]
                   - f_2 * osp1_222[k]
                   + f_3 * pc_x[k] * osd_444[k];

        t_741[k] = f_7 * osp0_223[k]
                   - f_8 * osp1_223[k]
                   + f_3 * pc_x[k] * osd_445[k];

        t_742[k] = f_7 * osp0_224[k]
                   - f_8 * osp1_224[k]
                   + f_3 * pc_x[k] * osd_446[k];

        t_743[k] = f_3 * pc_x[k] * osd_447[k];
    }

#pragma omp simd aligned(t_744, t_745, t_746, t_747, t_748, pc_x, pc_y, pc_z, nsd_375, \
                         nsd_381, nsd_383, osp0_223, osp1_223, osd_447, osd_448, \
                         osd_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_744[k] = f_3 * pc_x[k] * osd_448[k];

        t_745[k] = f_3 * pc_x[k] * osd_449[k];

        t_746[k] = f_12 * nsd_381[k]
                   + f_1 * osp0_223[k]
                   - f_2 * osp1_223[k]
                   + f_3 * pc_y[k] * osd_447[k];

        t_747[k] = f_11 * nsd_375[k]
                   + f_3 * pc_z[k] * osd_447[k];

        t_748[k] = f_12 * nsd_383[k]
                   + f_3 * pc_y[k] * osd_449[k];
    }

#pragma omp simd aligned(t_749, t_750, t_751, pc_x, pc_z, nsd_377, osp0_224, osp0_225, \
                         osp0_226, osp1_224, osp1_225, osp1_226, osd_449, osd_450, \
                         osd_451 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_749[k] = f_11 * nsd_377[k]
                   + f_1 * osp0_224[k]
                   - f_2 * osp1_224[k]
                   + f_3 * pc_z[k] * osd_449[k];

        t_750[k] = f_1 * osp0_225[k]
                   - f_2 * osp1_225[k]
                   + f_3 * pc_x[k] * osd_450[k];

        t_751[k] = f_7 * osp0_226[k]
                   - f_8 * osp1_226[k]
                   + f_3 * pc_x[k] * osd_451[k];
    }

#pragma omp simd aligned(t_752, t_753, t_754, t_755, t_756, pc_x, pc_y, nsd_387, osp0_226, \
                         osp0_227, osp1_226, osp1_227, osd_452, osd_453, osd_454, \
                         osd_455 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_752[k] = f_7 * osp0_227[k]
                   - f_8 * osp1_227[k]
                   + f_3 * pc_x[k] * osd_452[k];

        t_753[k] = f_3 * pc_x[k] * osd_453[k];

        t_754[k] = f_3 * pc_x[k] * osd_454[k];

        t_755[k] = f_3 * pc_x[k] * osd_455[k];

        t_756[k] = f_10 * nsd_387[k]
                   + f_1 * osp0_226[k]
                   - f_2 * osp1_226[k]
                   + f_3 * pc_y[k] * osd_453[k];
    }
}

static auto
compute_prim_osf_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsf0,
                                                          const size_t nsd, const size_t nsf1,
                                                          const size_t osp0, const size_t osp1,
                                                          const size_t osd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 5.0 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_9 = 4.5 / q;
    const auto f_10 = 1.0 / q;
    const auto f_12 = 1.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsf0_650 = buffer.data(nsf0 + 650);
    const auto *nsf0_652 = buffer.data(nsf0 + 652);
    const auto *nsf0_656 = buffer.data(nsf0 + 656);
    const auto *nsf0_659 = buffer.data(nsf0 + 659);

    const auto *nsd_381 = buffer.data(nsd + 381);
    const auto *nsd_383 = buffer.data(nsd + 383);
    const auto *nsd_387 = buffer.data(nsd + 387);
    const auto *nsd_389 = buffer.data(nsd + 389);
    const auto *nsd_393 = buffer.data(nsd + 393);
    const auto *nsd_395 = buffer.data(nsd + 395);

    const auto *nsf1_650 = buffer.data(nsf1 + 650);
    const auto *nsf1_652 = buffer.data(nsf1 + 652);
    const auto *nsf1_656 = buffer.data(nsf1 + 656);
    const auto *nsf1_659 = buffer.data(nsf1 + 659);

    const auto *osp0_227 = buffer.data(osp0 + 227);
    const auto *osp0_229 = buffer.data(osp0 + 229);
    const auto *osp0_231 = buffer.data(osp0 + 231);
    const auto *osp0_232 = buffer.data(osp0 + 232);
    const auto *osp0_233 = buffer.data(osp0 + 233);

    const auto *osp1_227 = buffer.data(osp1 + 227);
    const auto *osp1_229 = buffer.data(osp1 + 229);
    const auto *osp1_231 = buffer.data(osp1 + 231);
    const auto *osp1_232 = buffer.data(osp1 + 232);
    const auto *osp1_233 = buffer.data(osp1 + 233);

    const auto *osd_453 = buffer.data(osd + 453);
    const auto *osd_455 = buffer.data(osd + 455);
    const auto *osd_457 = buffer.data(osd + 457);
    const auto *osd_459 = buffer.data(osd + 459);
    const auto *osd_460 = buffer.data(osd + 460);
    const auto *osd_461 = buffer.data(osd + 461);
    const auto *osd_462 = buffer.data(osd + 462);
    const auto *osd_464 = buffer.data(osd + 464);
    const auto *osd_465 = buffer.data(osd + 465);
    const auto *osd_466 = buffer.data(osd + 466);
    const auto *osd_467 = buffer.data(osd + 467);

#pragma omp simd aligned(t_757, t_758, t_759, t_760, pa_y, pc_y, pc_z, nsf0_650, nsd_381, \
                         nsd_383, nsd_389, nsf1_650, osp0_227, osp1_227, osd_453, \
                         osd_455 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_757[k] = f_9 * nsd_381[k]
                   + f_3 * pc_z[k] * osd_453[k];

        t_758[k] = f_10 * nsd_389[k]
                   + f_3 * pc_y[k] * osd_455[k];

        t_759[k] = f_9 * nsd_383[k]
                   + f_1 * osp0_227[k]
                   - f_2 * osp1_227[k]
                   + f_3 * pc_z[k] * osd_455[k];

        t_760[k] = pa_y[k] * nsf0_650[k]
                   - f_4 * pc_y[k] * nsf1_650[k];
    }

#pragma omp simd aligned(t_761, t_762, t_763, t_764, t_765, pa_y, pc_x, pc_y, nsf0_652, \
                         nsf1_652, osp0_229, osp1_229, osd_457, osd_459, osd_460, \
                         osd_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_761[k] = f_7 * osp0_229[k]
                   - f_8 * osp1_229[k]
                   + f_3 * pc_x[k] * osd_457[k];

        t_762[k] = pa_y[k] * nsf0_652[k]
                   - f_4 * pc_y[k] * nsf1_652[k];

        t_763[k] = f_3 * pc_x[k] * osd_459[k];

        t_764[k] = f_3 * pc_x[k] * osd_460[k];

        t_765[k] = f_3 * pc_x[k] * osd_461[k];
    }

#pragma omp simd aligned(t_766, t_767, t_768, t_769, pa_y, pc_y, pc_z, nsf0_656, nsf0_659, \
                         nsd_387, nsd_393, nsd_395, nsf1_656, nsf1_659, osd_459, \
                         osd_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_766[k] = pa_y[k] * nsf0_656[k]
                   + f_12 * nsd_393[k]
                   - f_4 * pc_y[k] * nsf1_656[k];

        t_767[k] = f_6 * nsd_387[k]
                   + f_3 * pc_z[k] * osd_459[k];

        t_768[k] = f_5 * nsd_395[k]
                   + f_3 * pc_y[k] * osd_461[k];

        t_769[k] = pa_y[k] * nsf0_659[k]
                   - f_4 * pc_y[k] * nsf1_659[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, t_773, t_774, pc_x, pc_y, osp0_231, osp0_233, \
                         osp1_231, osp1_233, osd_462, osd_464, osd_465, \
                         osd_466 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = f_1 * osp0_231[k]
                   - f_2 * osp1_231[k]
                   + f_3 * pc_x[k] * osd_462[k];

        t_771[k] = f_3 * pc_y[k] * osd_462[k];

        t_772[k] = f_7 * osp0_233[k]
                   - f_8 * osp1_233[k]
                   + f_3 * pc_x[k] * osd_464[k];

        t_773[k] = f_3 * pc_x[k] * osd_465[k];

        t_774[k] = f_3 * pc_x[k] * osd_466[k];
    }

#pragma omp simd aligned(t_775, t_776, t_777, t_778, t_779, pc_x, pc_y, pc_z, nsd_395, \
                         osp0_232, osp0_233, osp1_232, osp1_233, osd_465, osd_466, \
                         osd_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_775[k] = f_3 * pc_x[k] * osd_467[k];

        t_776[k] = f_1 * osp0_232[k]
                   - f_2 * osp1_232[k]
                   + f_3 * pc_y[k] * osd_465[k];

        t_777[k] = f_7 * osp0_233[k]
                   - f_8 * osp1_233[k]
                   + f_3 * pc_y[k] * osd_466[k];

        t_778[k] = f_3 * pc_y[k] * osd_467[k];

        t_779[k] = f_0 * nsd_395[k]
                   + f_1 * osp0_233[k]
                   - f_2 * osp1_233[k]
                   + f_3 * pc_z[k] * osd_467[k];
    }
}

auto
compute_prim_osf_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t nsf0, const size_t nsd,
                                                   const size_t nsf1, const size_t osp0,
                                                   const size_t osp1, const size_t osd,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_osf_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, nsf0, nsd,
                                                              nsf1, osp0, osp1, osd, ncols,
                                                              gamma, p, q);

    compute_prim_osf_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, nsf0, nsd,
                                                              nsf1, osp0, osp1, osd, ncols,
                                                              gamma, p, q);

    compute_prim_osf_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, nsf0, nsd,
                                                              nsf1, osp0, osp1, osd, ncols,
                                                              gamma, p, q);

    compute_prim_osf_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, nsf0, nsd,
                                                              nsf1, osp0, osp1, osd, ncols,
                                                              gamma, p, q);

    compute_prim_osf_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, nsf0, nsd,
                                                              nsf1, osp0, osp1, osd, ncols,
                                                              gamma, p, q);

    compute_prim_osf_three_center_electron_repulsion_0_piece5(buffer, target, pa, pc, nsf0, nsd,
                                                              nsf1, osp0, osp1, osd, ncols,
                                                              gamma, p, q);

    compute_prim_osf_three_center_electron_repulsion_0_piece6(buffer, target, pa, pc, nsf0, nsd,
                                                              nsf1, osp0, osp1, osd, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
