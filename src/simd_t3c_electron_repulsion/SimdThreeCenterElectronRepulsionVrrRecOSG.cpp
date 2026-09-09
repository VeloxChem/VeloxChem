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


#include "SimdThreeCenterElectronRepulsionVrrRecOSG.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_osg_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsg0,
                                                          const size_t nsf, const size_t nsg1,
                                                          const size_t osd0, const size_t osd1,
                                                          const size_t osf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 5.0 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 4.5 / q;
    const auto f_13 = 4.0 / q;
    const auto f_14 = 1.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsg0_0 = buffer.data(nsg0 + 0);
    const auto *nsg0_3 = buffer.data(nsg0 + 3);
    const auto *nsg0_5 = buffer.data(nsg0 + 5);
    const auto *nsg0_10 = buffer.data(nsg0 + 10);
    const auto *nsg0_14 = buffer.data(nsg0 + 14);
    const auto *nsg0_18 = buffer.data(nsg0 + 18);
    const auto *nsg0_25 = buffer.data(nsg0 + 25);
    const auto *nsg0_30 = buffer.data(nsg0 + 30);
    const auto *nsg0_35 = buffer.data(nsg0 + 35);
    const auto *nsg0_44 = buffer.data(nsg0 + 44);
    const auto *nsg0_45 = buffer.data(nsg0 + 45);
    const auto *nsg0_48 = buffer.data(nsg0 + 48);
    const auto *nsg0_55 = buffer.data(nsg0 + 55);
    const auto *nsg0_75 = buffer.data(nsg0 + 75);
    const auto *nsg0_78 = buffer.data(nsg0 + 78);
    const auto *nsg0_80 = buffer.data(nsg0 + 80);

    const auto *nsf_0 = buffer.data(nsf + 0);
    const auto *nsf_1 = buffer.data(nsf + 1);
    const auto *nsf_2 = buffer.data(nsf + 2);
    const auto *nsf_6 = buffer.data(nsf + 6);
    const auto *nsf_9 = buffer.data(nsf + 9);
    const auto *nsf_10 = buffer.data(nsf + 10);
    const auto *nsf_16 = buffer.data(nsf + 16);
    const auto *nsf_18 = buffer.data(nsf + 18);
    const auto *nsf_19 = buffer.data(nsf + 19);
    const auto *nsf_20 = buffer.data(nsf + 20);
    const auto *nsf_22 = buffer.data(nsf + 22);
    const auto *nsf_26 = buffer.data(nsf + 26);
    const auto *nsf_27 = buffer.data(nsf + 27);
    const auto *nsf_28 = buffer.data(nsf + 28);
    const auto *nsf_29 = buffer.data(nsf + 29);
    const auto *nsf_30 = buffer.data(nsf + 30);
    const auto *nsf_33 = buffer.data(nsf + 33);
    const auto *nsf_36 = buffer.data(nsf + 36);
    const auto *nsf_38 = buffer.data(nsf + 38);
    const auto *nsf_39 = buffer.data(nsf + 39);
    const auto *nsf_40 = buffer.data(nsf + 40);
    const auto *nsf_42 = buffer.data(nsf + 42);
    const auto *nsf_46 = buffer.data(nsf + 46);
    const auto *nsf_47 = buffer.data(nsf + 47);
    const auto *nsf_48 = buffer.data(nsf + 48);
    const auto *nsf_49 = buffer.data(nsf + 49);
    const auto *nsf_50 = buffer.data(nsf + 50);
    const auto *nsf_51 = buffer.data(nsf + 51);
    const auto *nsf_52 = buffer.data(nsf + 52);
    const auto *nsf_55 = buffer.data(nsf + 55);
    const auto *nsf_56 = buffer.data(nsf + 56);
    const auto *nsf_57 = buffer.data(nsf + 57);
    const auto *nsf_59 = buffer.data(nsf + 59);
    const auto *nsf_60 = buffer.data(nsf + 60);
    const auto *nsf_63 = buffer.data(nsf + 63);
    const auto *nsf_66 = buffer.data(nsf + 66);
    const auto *nsf_68 = buffer.data(nsf + 68);
    const auto *nsf_69 = buffer.data(nsf + 69);
    const auto *nsf_75 = buffer.data(nsf + 75);
    const auto *nsf_76 = buffer.data(nsf + 76);
    const auto *nsf_77 = buffer.data(nsf + 77);
    const auto *nsf_78 = buffer.data(nsf + 78);
    const auto *nsf_79 = buffer.data(nsf + 79);
    const auto *nsf_86 = buffer.data(nsf + 86);
    const auto *nsf_87 = buffer.data(nsf + 87);
    const auto *nsf_88 = buffer.data(nsf + 88);
    const auto *nsf_89 = buffer.data(nsf + 89);

    const auto *nsg1_0 = buffer.data(nsg1 + 0);
    const auto *nsg1_3 = buffer.data(nsg1 + 3);
    const auto *nsg1_5 = buffer.data(nsg1 + 5);
    const auto *nsg1_10 = buffer.data(nsg1 + 10);
    const auto *nsg1_14 = buffer.data(nsg1 + 14);
    const auto *nsg1_18 = buffer.data(nsg1 + 18);
    const auto *nsg1_25 = buffer.data(nsg1 + 25);
    const auto *nsg1_30 = buffer.data(nsg1 + 30);
    const auto *nsg1_35 = buffer.data(nsg1 + 35);
    const auto *nsg1_44 = buffer.data(nsg1 + 44);
    const auto *nsg1_45 = buffer.data(nsg1 + 45);
    const auto *nsg1_48 = buffer.data(nsg1 + 48);
    const auto *nsg1_55 = buffer.data(nsg1 + 55);
    const auto *nsg1_75 = buffer.data(nsg1 + 75);
    const auto *nsg1_78 = buffer.data(nsg1 + 78);
    const auto *nsg1_80 = buffer.data(nsg1 + 80);

    const auto *osd0_0 = buffer.data(osd0 + 0);
    const auto *osd0_3 = buffer.data(osd0 + 3);
    const auto *osd0_5 = buffer.data(osd0 + 5);
    const auto *osd0_9 = buffer.data(osd0 + 9);
    const auto *osd0_16 = buffer.data(osd0 + 16);
    const auto *osd0_17 = buffer.data(osd0 + 17);
    const auto *osd0_18 = buffer.data(osd0 + 18);
    const auto *osd0_21 = buffer.data(osd0 + 21);
    const auto *osd0_23 = buffer.data(osd0 + 23);
    const auto *osd0_29 = buffer.data(osd0 + 29);
    const auto *osd0_30 = buffer.data(osd0 + 30);
    const auto *osd0_33 = buffer.data(osd0 + 33);
    const auto *osd0_34 = buffer.data(osd0 + 34);
    const auto *osd0_35 = buffer.data(osd0 + 35);
    const auto *osd0_36 = buffer.data(osd0 + 36);
    const auto *osd0_39 = buffer.data(osd0 + 39);
    const auto *osd0_41 = buffer.data(osd0 + 41);
    const auto *osd0_47 = buffer.data(osd0 + 47);

    const auto *osd1_0 = buffer.data(osd1 + 0);
    const auto *osd1_3 = buffer.data(osd1 + 3);
    const auto *osd1_5 = buffer.data(osd1 + 5);
    const auto *osd1_9 = buffer.data(osd1 + 9);
    const auto *osd1_16 = buffer.data(osd1 + 16);
    const auto *osd1_17 = buffer.data(osd1 + 17);
    const auto *osd1_18 = buffer.data(osd1 + 18);
    const auto *osd1_21 = buffer.data(osd1 + 21);
    const auto *osd1_23 = buffer.data(osd1 + 23);
    const auto *osd1_29 = buffer.data(osd1 + 29);
    const auto *osd1_30 = buffer.data(osd1 + 30);
    const auto *osd1_33 = buffer.data(osd1 + 33);
    const auto *osd1_34 = buffer.data(osd1 + 34);
    const auto *osd1_35 = buffer.data(osd1 + 35);
    const auto *osd1_36 = buffer.data(osd1 + 36);
    const auto *osd1_39 = buffer.data(osd1 + 39);
    const auto *osd1_41 = buffer.data(osd1 + 41);
    const auto *osd1_47 = buffer.data(osd1 + 47);

    const auto *osf_0 = buffer.data(osf + 0);
    const auto *osf_1 = buffer.data(osf + 1);
    const auto *osf_2 = buffer.data(osf + 2);
    const auto *osf_3 = buffer.data(osf + 3);
    const auto *osf_5 = buffer.data(osf + 5);
    const auto *osf_6 = buffer.data(osf + 6);
    const auto *osf_8 = buffer.data(osf + 8);
    const auto *osf_9 = buffer.data(osf + 9);
    const auto *osf_10 = buffer.data(osf + 10);
    const auto *osf_11 = buffer.data(osf + 11);
    const auto *osf_13 = buffer.data(osf + 13);
    const auto *osf_16 = buffer.data(osf + 16);
    const auto *osf_17 = buffer.data(osf + 17);
    const auto *osf_18 = buffer.data(osf + 18);
    const auto *osf_19 = buffer.data(osf + 19);
    const auto *osf_20 = buffer.data(osf + 20);
    const auto *osf_22 = buffer.data(osf + 22);
    const auto *osf_25 = buffer.data(osf + 25);
    const auto *osf_26 = buffer.data(osf + 26);
    const auto *osf_27 = buffer.data(osf + 27);
    const auto *osf_28 = buffer.data(osf + 28);
    const auto *osf_29 = buffer.data(osf + 29);
    const auto *osf_30 = buffer.data(osf + 30);
    const auto *osf_31 = buffer.data(osf + 31);
    const auto *osf_32 = buffer.data(osf + 32);
    const auto *osf_33 = buffer.data(osf + 33);
    const auto *osf_36 = buffer.data(osf + 36);
    const auto *osf_37 = buffer.data(osf + 37);
    const auto *osf_38 = buffer.data(osf + 38);
    const auto *osf_39 = buffer.data(osf + 39);
    const auto *osf_40 = buffer.data(osf + 40);
    const auto *osf_42 = buffer.data(osf + 42);
    const auto *osf_46 = buffer.data(osf + 46);
    const auto *osf_47 = buffer.data(osf + 47);
    const auto *osf_48 = buffer.data(osf + 48);
    const auto *osf_49 = buffer.data(osf + 49);
    const auto *osf_50 = buffer.data(osf + 50);
    const auto *osf_51 = buffer.data(osf + 51);
    const auto *osf_52 = buffer.data(osf + 52);
    const auto *osf_55 = buffer.data(osf + 55);
    const auto *osf_56 = buffer.data(osf + 56);
    const auto *osf_57 = buffer.data(osf + 57);
    const auto *osf_58 = buffer.data(osf + 58);
    const auto *osf_59 = buffer.data(osf + 59);
    const auto *osf_60 = buffer.data(osf + 60);
    const auto *osf_61 = buffer.data(osf + 61);
    const auto *osf_62 = buffer.data(osf + 62);
    const auto *osf_63 = buffer.data(osf + 63);
    const auto *osf_66 = buffer.data(osf + 66);
    const auto *osf_67 = buffer.data(osf + 67);
    const auto *osf_68 = buffer.data(osf + 68);
    const auto *osf_69 = buffer.data(osf + 69);
    const auto *osf_70 = buffer.data(osf + 70);
    const auto *osf_72 = buffer.data(osf + 72);
    const auto *osf_75 = buffer.data(osf + 75);
    const auto *osf_76 = buffer.data(osf + 76);
    const auto *osf_77 = buffer.data(osf + 77);
    const auto *osf_78 = buffer.data(osf + 78);
    const auto *osf_79 = buffer.data(osf + 79);
    const auto *osf_80 = buffer.data(osf + 80);
    const auto *osf_82 = buffer.data(osf + 82);
    const auto *osf_86 = buffer.data(osf + 86);
    const auto *osf_87 = buffer.data(osf + 87);
    const auto *osf_88 = buffer.data(osf + 88);
    const auto *osf_89 = buffer.data(osf + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, nsf_0, osd0_0, \
                         osd1_0, osf_0, osf_1, osf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * nsf_0[k]
                 + f_1 * osd0_0[k]
                 - f_2 * osd1_0[k]
                 + f_3 * pc_x[k] * osf_0[k];

        t_1[k] = f_3 * pc_y[k] * osf_0[k];

        t_2[k] = f_3 * pc_z[k] * osf_0[k];

        t_3[k] = f_4 * osd0_0[k]
                 - f_5 * osd1_0[k]
                 + f_3 * pc_y[k] * osf_1[k];

        t_4[k] = f_3 * pc_y[k] * osf_2[k];

        t_5[k] = f_4 * osd0_0[k]
                 - f_5 * osd1_0[k]
                 + f_3 * pc_z[k] * osf_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_x, pc_y, pc_z, nsf_6, nsf_9, osd0_3, \
                         osd1_3, osf_3, osf_5, osf_6, osf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * nsf_6[k]
                 + f_3 * pc_x[k] * osf_6[k];

        t_7[k] = f_3 * pc_z[k] * osf_3[k];

        t_8[k] = f_3 * pc_y[k] * osf_5[k];

        t_9[k] = f_0 * nsf_9[k]
                 + f_3 * pc_x[k] * osf_9[k];

        t_10[k] = f_1 * osd0_3[k]
                  - f_2 * osd1_3[k]
                  + f_3 * pc_y[k] * osf_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pa_y, pc_y, pc_z, nsg0_0, nsg1_0, \
                         osd0_5, osd1_5, osf_6, osf_8, osf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * osf_6[k];

        t_12[k] = f_4 * osd0_5[k]
                  - f_5 * osd1_5[k]
                  + f_3 * pc_y[k] * osf_8[k];

        t_13[k] = f_3 * pc_y[k] * osf_9[k];

        t_14[k] = f_1 * osd0_5[k]
                  - f_2 * osd1_5[k]
                  + f_3 * pc_z[k] * osf_9[k];

        t_15[k] = pa_y[k] * nsg0_0[k]
                  - f_6 * pc_y[k] * nsg1_0[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pa_y, pc_y, pc_z, nsg0_3, nsg0_5, \
                         nsf_0, nsf_1, nsg1_3, nsg1_5, osf_10, osf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_7 * nsf_0[k]
                  + f_3 * pc_y[k] * osf_10[k];

        t_17[k] = f_3 * pc_z[k] * osf_10[k];

        t_18[k] = pa_y[k] * nsg0_3[k]
                  + f_8 * nsf_1[k]
                  - f_6 * pc_y[k] * nsg1_3[k];

        t_19[k] = f_3 * pc_z[k] * osf_11[k];

        t_20[k] = pa_y[k] * nsg0_5[k]
                  - f_6 * pc_y[k] * nsg1_5[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_x, pc_z, nsf_16, nsf_18, nsf_19, osf_13, \
                         osf_16, osf_18, osf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_9 * nsf_16[k]
                  + f_3 * pc_x[k] * osf_16[k];

        t_22[k] = f_3 * pc_z[k] * osf_13[k];

        t_23[k] = f_9 * nsf_18[k]
                  + f_3 * pc_x[k] * osf_18[k];

        t_24[k] = f_9 * nsf_19[k]
                  + f_3 * pc_x[k] * osf_19[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pc_y, pc_z, nsf_6, nsf_9, osd0_9, osd1_9, \
                         osf_16, osf_17, osf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_7 * nsf_6[k]
                  + f_1 * osd0_9[k]
                  - f_2 * osd1_9[k]
                  + f_3 * pc_y[k] * osf_16[k];

        t_26[k] = f_3 * pc_z[k] * osf_16[k];

        t_27[k] = f_4 * osd0_9[k]
                  - f_5 * osd1_9[k]
                  + f_3 * pc_z[k] * osf_17[k];

        t_28[k] = f_7 * nsf_9[k]
                  + f_3 * pc_y[k] * osf_19[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_y, pa_z, pc_y, pc_z, nsg0_0, nsg0_14, \
                         nsf_0, nsg1_0, nsg1_14, osf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_y[k] * nsg0_14[k]
                  - f_6 * pc_y[k] * nsg1_14[k];

        t_30[k] = pa_z[k] * nsg0_0[k]
                  - f_6 * pc_z[k] * nsg1_0[k];

        t_31[k] = f_3 * pc_y[k] * osf_20[k];

        t_32[k] = f_7 * nsf_0[k]
                  + f_3 * pc_z[k] * osf_20[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_z, pc_x, pc_y, pc_z, nsg0_3, nsg0_5, \
                         nsf_2, nsf_26, nsg1_3, nsg1_5, osf_22, \
                         osf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_z[k] * nsg0_3[k]
                  - f_6 * pc_z[k] * nsg1_3[k];

        t_34[k] = f_3 * pc_y[k] * osf_22[k];

        t_35[k] = pa_z[k] * nsg0_5[k]
                  + f_8 * nsf_2[k]
                  - f_6 * pc_z[k] * nsg1_5[k];

        t_36[k] = f_9 * nsf_26[k]
                  + f_3 * pc_x[k] * osf_26[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_z, pc_x, pc_y, pc_z, nsg0_10, nsf_27, \
                         nsf_29, nsg1_10, osf_25, osf_27, osf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_9 * nsf_27[k]
                  + f_3 * pc_x[k] * osf_27[k];

        t_38[k] = f_3 * pc_y[k] * osf_25[k];

        t_39[k] = f_9 * nsf_29[k]
                  + f_3 * pc_x[k] * osf_29[k];

        t_40[k] = pa_z[k] * nsg0_10[k]
                  - f_6 * pc_z[k] * nsg1_10[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pc_y, pc_z, nsf_9, osd0_16, osd0_17, osd1_16, \
                         osd1_17, osf_27, osf_28, osf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_10 * osd0_16[k]
                  - f_11 * osd1_16[k]
                  + f_3 * pc_y[k] * osf_27[k];

        t_42[k] = f_4 * osd0_17[k]
                  - f_5 * osd1_17[k]
                  + f_3 * pc_y[k] * osf_28[k];

        t_43[k] = f_3 * pc_y[k] * osf_29[k];

        t_44[k] = f_7 * nsf_9[k]
                  + f_1 * osd0_17[k]
                  - f_2 * osd1_17[k]
                  + f_3 * pc_z[k] * osf_29[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pc_x, pc_y, pc_z, nsf_10, nsf_30, nsf_33, \
                         osd0_18, osd0_21, osd1_18, osd1_21, osf_30, \
                         osf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_12 * nsf_30[k]
                  + f_1 * osd0_18[k]
                  - f_2 * osd1_18[k]
                  + f_3 * pc_x[k] * osf_30[k];

        t_46[k] = f_8 * nsf_10[k]
                  + f_3 * pc_y[k] * osf_30[k];

        t_47[k] = f_3 * pc_z[k] * osf_30[k];

        t_48[k] = f_12 * nsf_33[k]
                  + f_4 * osd0_21[k]
                  - f_5 * osd1_21[k]
                  + f_3 * pc_x[k] * osf_33[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, pc_x, pc_z, nsf_36, nsf_38, osd0_18, \
                         osd1_18, osf_31, osf_32, osf_33, osf_36, \
                         osf_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_3 * pc_z[k] * osf_31[k];

        t_50[k] = f_4 * osd0_18[k]
                  - f_5 * osd1_18[k]
                  + f_3 * pc_z[k] * osf_32[k];

        t_51[k] = f_12 * nsf_36[k]
                  + f_3 * pc_x[k] * osf_36[k];

        t_52[k] = f_3 * pc_z[k] * osf_33[k];

        t_53[k] = f_12 * nsf_38[k]
                  + f_3 * pc_x[k] * osf_38[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pc_x, pc_y, pc_z, nsf_16, nsf_19, \
                         nsf_39, osd0_21, osd1_21, osf_36, osf_37, \
                         osf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_12 * nsf_39[k]
                  + f_3 * pc_x[k] * osf_39[k];

        t_55[k] = f_8 * nsf_16[k]
                  + f_1 * osd0_21[k]
                  - f_2 * osd1_21[k]
                  + f_3 * pc_y[k] * osf_36[k];

        t_56[k] = f_3 * pc_z[k] * osf_36[k];

        t_57[k] = f_4 * osd0_21[k]
                  - f_5 * osd1_21[k]
                  + f_3 * pc_z[k] * osf_37[k];

        t_58[k] = f_8 * nsf_19[k]
                  + f_3 * pc_y[k] * osf_39[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_y, pc_y, pc_z, nsg0_30, nsf_10, nsf_20, \
                         nsg1_30, osd0_23, osd1_23, osf_39, osf_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_1 * osd0_23[k]
                  - f_2 * osd1_23[k]
                  + f_3 * pc_z[k] * osf_39[k];

        t_60[k] = pa_y[k] * nsg0_30[k]
                  - f_6 * pc_y[k] * nsg1_30[k];

        t_61[k] = f_7 * nsf_20[k]
                  + f_3 * pc_y[k] * osf_40[k];

        t_62[k] = f_7 * nsf_10[k]
                  + f_3 * pc_z[k] * osf_40[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, pa_y, pa_z, pc_y, pc_z, nsg0_18, nsg0_35, nsf_22, \
                         nsg1_18, nsg1_35, osf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = pa_z[k] * nsg0_18[k]
                  - f_6 * pc_z[k] * nsg1_18[k];

        t_64[k] = f_7 * nsf_22[k]
                  + f_3 * pc_y[k] * osf_42[k];

        t_65[k] = pa_y[k] * nsg0_35[k]
                  - f_6 * pc_y[k] * nsg1_35[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pc_x, nsf_46, nsf_47, nsf_48, nsf_49, osf_46, \
                         osf_47, osf_48, osf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_12 * nsf_46[k]
                  + f_3 * pc_x[k] * osf_46[k];

        t_67[k] = f_12 * nsf_47[k]
                  + f_3 * pc_x[k] * osf_47[k];

        t_68[k] = f_12 * nsf_48[k]
                  + f_3 * pc_x[k] * osf_48[k];

        t_69[k] = f_12 * nsf_49[k]
                  + f_3 * pc_x[k] * osf_49[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pa_z, pc_y, pc_z, nsg0_25, nsf_16, nsf_28, nsg1_25, \
                         osd0_29, osd1_29, osf_46, osf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pa_z[k] * nsg0_25[k]
                  - f_6 * pc_z[k] * nsg1_25[k];

        t_71[k] = f_7 * nsf_16[k]
                  + f_3 * pc_z[k] * osf_46[k];

        t_72[k] = f_7 * nsf_28[k]
                  + f_4 * osd0_29[k]
                  - f_5 * osd1_29[k]
                  + f_3 * pc_y[k] * osf_48[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_y, pc_x, pc_y, nsg0_44, nsf_29, nsf_50, \
                         nsg1_44, osd0_30, osd1_30, osf_49, osf_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_7 * nsf_29[k]
                  + f_3 * pc_y[k] * osf_49[k];

        t_74[k] = pa_y[k] * nsg0_44[k]
                  - f_6 * pc_y[k] * nsg1_44[k];

        t_75[k] = f_12 * nsf_50[k]
                  + f_1 * osd0_30[k]
                  - f_2 * osd1_30[k]
                  + f_3 * pc_x[k] * osf_50[k];

        t_76[k] = f_3 * pc_y[k] * osf_50[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pc_y, pc_z, nsf_20, osd0_30, osd1_30, osf_50, \
                         osf_51, osf_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_8 * nsf_20[k]
                  + f_3 * pc_z[k] * osf_50[k];

        t_78[k] = f_4 * osd0_30[k]
                  - f_5 * osd1_30[k]
                  + f_3 * pc_y[k] * osf_51[k];

        t_79[k] = f_3 * pc_y[k] * osf_52[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pc_x, pc_y, nsf_55, nsf_56, nsf_57, osd0_35, \
                         osd1_35, osf_55, osf_56, osf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_12 * nsf_55[k]
                  + f_4 * osd0_35[k]
                  - f_5 * osd1_35[k]
                  + f_3 * pc_x[k] * osf_55[k];

        t_81[k] = f_12 * nsf_56[k]
                  + f_3 * pc_x[k] * osf_56[k];

        t_82[k] = f_12 * nsf_57[k]
                  + f_3 * pc_x[k] * osf_57[k];

        t_83[k] = f_3 * pc_y[k] * osf_55[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pc_x, pc_y, nsf_59, osd0_33, osd0_34, osd1_33, \
                         osd1_34, osf_56, osf_57, osf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_12 * nsf_59[k]
                  + f_3 * pc_x[k] * osf_59[k];

        t_85[k] = f_1 * osd0_33[k]
                  - f_2 * osd1_33[k]
                  + f_3 * pc_y[k] * osf_56[k];

        t_86[k] = f_10 * osd0_34[k]
                  - f_11 * osd1_34[k]
                  + f_3 * pc_y[k] * osf_57[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pc_x, pc_y, pc_z, nsf_29, nsf_60, osd0_35, \
                         osd0_36, osd1_35, osd1_36, osf_58, osf_59, \
                         osf_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_4 * osd0_35[k]
                  - f_5 * osd1_35[k]
                  + f_3 * pc_y[k] * osf_58[k];

        t_88[k] = f_3 * pc_y[k] * osf_59[k];

        t_89[k] = f_8 * nsf_29[k]
                  + f_1 * osd0_35[k]
                  - f_2 * osd1_35[k]
                  + f_3 * pc_z[k] * osf_59[k];

        t_90[k] = f_13 * nsf_60[k]
                  + f_1 * osd0_36[k]
                  - f_2 * osd1_36[k]
                  + f_3 * pc_x[k] * osf_60[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pc_x, pc_y, pc_z, nsf_30, nsf_63, osd0_39, \
                         osd1_39, osf_60, osf_61, osf_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_14 * nsf_30[k]
                  + f_3 * pc_y[k] * osf_60[k];

        t_92[k] = f_3 * pc_z[k] * osf_60[k];

        t_93[k] = f_13 * nsf_63[k]
                  + f_4 * osd0_39[k]
                  - f_5 * osd1_39[k]
                  + f_3 * pc_x[k] * osf_63[k];

        t_94[k] = f_3 * pc_z[k] * osf_61[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pc_x, pc_z, nsf_66, nsf_68, osd0_36, osd1_36, \
                         osf_62, osf_63, osf_66, osf_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_4 * osd0_36[k]
                  - f_5 * osd1_36[k]
                  + f_3 * pc_z[k] * osf_62[k];

        t_96[k] = f_13 * nsf_66[k]
                  + f_3 * pc_x[k] * osf_66[k];

        t_97[k] = f_3 * pc_z[k] * osf_63[k];

        t_98[k] = f_13 * nsf_68[k]
                  + f_3 * pc_x[k] * osf_68[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, pc_x, pc_y, pc_z, nsf_36, nsf_39, \
                         nsf_69, osd0_39, osd1_39, osf_66, osf_67, \
                         osf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_13 * nsf_69[k]
                  + f_3 * pc_x[k] * osf_69[k];

        t_100[k] = f_14 * nsf_36[k]
                   + f_1 * osd0_39[k]
                   - f_2 * osd1_39[k]
                   + f_3 * pc_y[k] * osf_66[k];

        t_101[k] = f_3 * pc_z[k] * osf_66[k];

        t_102[k] = f_4 * osd0_39[k]
                   - f_5 * osd1_39[k]
                   + f_3 * pc_z[k] * osf_67[k];

        t_103[k] = f_14 * nsf_39[k]
                   + f_3 * pc_y[k] * osf_69[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pa_z, pc_y, pc_z, nsg0_45, nsf_30, \
                         nsf_40, nsg1_45, osd0_41, osd1_41, osf_69, \
                         osf_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_1 * osd0_41[k]
                   - f_2 * osd1_41[k]
                   + f_3 * pc_z[k] * osf_69[k];

        t_105[k] = pa_z[k] * nsg0_45[k]
                   - f_6 * pc_z[k] * nsg1_45[k];

        t_106[k] = f_8 * nsf_40[k]
                   + f_3 * pc_y[k] * osf_70[k];

        t_107[k] = f_7 * nsf_30[k]
                   + f_3 * pc_z[k] * osf_70[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pa_z, pc_x, pc_y, pc_z, nsg0_48, nsf_42, nsf_75, \
                         nsg1_48, osd0_47, osd1_47, osf_72, osf_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = pa_z[k] * nsg0_48[k]
                   - f_6 * pc_z[k] * nsg1_48[k];

        t_109[k] = f_8 * nsf_42[k]
                   + f_3 * pc_y[k] * osf_72[k];

        t_110[k] = f_13 * nsf_75[k]
                   + f_4 * osd0_47[k]
                   - f_5 * osd1_47[k]
                   + f_3 * pc_x[k] * osf_75[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pc_x, nsf_76, nsf_77, nsf_78, nsf_79, \
                         osf_76, osf_77, osf_78, osf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_13 * nsf_76[k]
                   + f_3 * pc_x[k] * osf_76[k];

        t_112[k] = f_13 * nsf_77[k]
                   + f_3 * pc_x[k] * osf_77[k];

        t_113[k] = f_13 * nsf_78[k]
                   + f_3 * pc_x[k] * osf_78[k];

        t_114[k] = f_13 * nsf_79[k]
                   + f_3 * pc_x[k] * osf_79[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, pa_z, pc_y, pc_z, nsg0_55, nsf_36, nsf_48, \
                         nsg1_55, osd0_47, osd1_47, osf_76, osf_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = pa_z[k] * nsg0_55[k]
                   - f_6 * pc_z[k] * nsg1_55[k];

        t_116[k] = f_7 * nsf_36[k]
                   + f_3 * pc_z[k] * osf_76[k];

        t_117[k] = f_8 * nsf_48[k]
                   + f_4 * osd0_47[k]
                   - f_5 * osd1_47[k]
                   + f_3 * pc_y[k] * osf_78[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, t_121, pa_y, pc_y, pc_z, nsg0_75, nsf_39, \
                         nsf_49, nsf_50, nsg1_75, osd0_47, osd1_47, osf_79, \
                         osf_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_8 * nsf_49[k]
                   + f_3 * pc_y[k] * osf_79[k];

        t_119[k] = f_7 * nsf_39[k]
                   + f_1 * osd0_47[k]
                   - f_2 * osd1_47[k]
                   + f_3 * pc_z[k] * osf_79[k];

        t_120[k] = pa_y[k] * nsg0_75[k]
                   - f_6 * pc_y[k] * nsg1_75[k];

        t_121[k] = f_7 * nsf_50[k]
                   + f_3 * pc_y[k] * osf_80[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, pa_y, pc_y, pc_z, nsg0_78, nsg0_80, \
                         nsf_40, nsf_51, nsf_52, nsg1_78, nsg1_80, osf_80, \
                         osf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_8 * nsf_40[k]
                   + f_3 * pc_z[k] * osf_80[k];

        t_123[k] = pa_y[k] * nsg0_78[k]
                   + f_8 * nsf_51[k]
                   - f_6 * pc_y[k] * nsg1_78[k];

        t_124[k] = f_7 * nsf_52[k]
                   + f_3 * pc_y[k] * osf_82[k];

        t_125[k] = pa_y[k] * nsg0_80[k]
                   - f_6 * pc_y[k] * nsg1_80[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pc_x, nsf_86, nsf_87, nsf_88, nsf_89, \
                         osf_86, osf_87, osf_88, osf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_13 * nsf_86[k]
                   + f_3 * pc_x[k] * osf_86[k];

        t_127[k] = f_13 * nsf_87[k]
                   + f_3 * pc_x[k] * osf_87[k];

        t_128[k] = f_13 * nsf_88[k]
                   + f_3 * pc_x[k] * osf_88[k];

        t_129[k] = f_13 * nsf_89[k]
                   + f_3 * pc_x[k] * osf_89[k];
    }
}

static auto
compute_prim_osg_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsg0,
                                                          const size_t nsf, const size_t nsg1,
                                                          const size_t osd0, const size_t osd1,
                                                          const size_t osf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_13 = 4.0 / q;
    const auto f_14 = 1.5 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.5 / q;

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
    auto *t_249 = buffer.data(target + 249);
    auto *t_250 = buffer.data(target + 250);
    auto *t_251 = buffer.data(target + 251);
    auto *t_252 = buffer.data(target + 252);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsg0_89 = buffer.data(nsg0 + 89);
    const auto *nsg0_90 = buffer.data(nsg0 + 90);
    const auto *nsg0_93 = buffer.data(nsg0 + 93);
    const auto *nsg0_100 = buffer.data(nsg0 + 100);
    const auto *nsg0_135 = buffer.data(nsg0 + 135);
    const auto *nsg0_138 = buffer.data(nsg0 + 138);
    const auto *nsg0_140 = buffer.data(nsg0 + 140);
    const auto *nsg0_149 = buffer.data(nsg0 + 149);
    const auto *nsg0_150 = buffer.data(nsg0 + 150);
    const auto *nsg0_153 = buffer.data(nsg0 + 153);
    const auto *nsg0_160 = buffer.data(nsg0 + 160);

    const auto *nsf_46 = buffer.data(nsf + 46);
    const auto *nsf_50 = buffer.data(nsf + 50);
    const auto *nsf_56 = buffer.data(nsf + 56);
    const auto *nsf_58 = buffer.data(nsf + 58);
    const auto *nsf_59 = buffer.data(nsf + 59);
    const auto *nsf_60 = buffer.data(nsf + 60);
    const auto *nsf_66 = buffer.data(nsf + 66);
    const auto *nsf_69 = buffer.data(nsf + 69);
    const auto *nsf_70 = buffer.data(nsf + 70);
    const auto *nsf_72 = buffer.data(nsf + 72);
    const auto *nsf_76 = buffer.data(nsf + 76);
    const auto *nsf_78 = buffer.data(nsf + 78);
    const auto *nsf_79 = buffer.data(nsf + 79);
    const auto *nsf_80 = buffer.data(nsf + 80);
    const auto *nsf_82 = buffer.data(nsf + 82);
    const auto *nsf_86 = buffer.data(nsf + 86);
    const auto *nsf_88 = buffer.data(nsf + 88);
    const auto *nsf_89 = buffer.data(nsf + 89);
    const auto *nsf_90 = buffer.data(nsf + 90);
    const auto *nsf_91 = buffer.data(nsf + 91);
    const auto *nsf_92 = buffer.data(nsf + 92);
    const auto *nsf_95 = buffer.data(nsf + 95);
    const auto *nsf_96 = buffer.data(nsf + 96);
    const auto *nsf_97 = buffer.data(nsf + 97);
    const auto *nsf_98 = buffer.data(nsf + 98);
    const auto *nsf_99 = buffer.data(nsf + 99);
    const auto *nsf_100 = buffer.data(nsf + 100);
    const auto *nsf_103 = buffer.data(nsf + 103);
    const auto *nsf_106 = buffer.data(nsf + 106);
    const auto *nsf_108 = buffer.data(nsf + 108);
    const auto *nsf_109 = buffer.data(nsf + 109);
    const auto *nsf_110 = buffer.data(nsf + 110);
    const auto *nsf_112 = buffer.data(nsf + 112);
    const auto *nsf_115 = buffer.data(nsf + 115);
    const auto *nsf_116 = buffer.data(nsf + 116);
    const auto *nsf_117 = buffer.data(nsf + 117);
    const auto *nsf_118 = buffer.data(nsf + 118);
    const auto *nsf_119 = buffer.data(nsf + 119);
    const auto *nsf_120 = buffer.data(nsf + 120);
    const auto *nsf_123 = buffer.data(nsf + 123);
    const auto *nsf_125 = buffer.data(nsf + 125);
    const auto *nsf_126 = buffer.data(nsf + 126);
    const auto *nsf_127 = buffer.data(nsf + 127);
    const auto *nsf_128 = buffer.data(nsf + 128);
    const auto *nsf_129 = buffer.data(nsf + 129);
    const auto *nsf_136 = buffer.data(nsf + 136);
    const auto *nsf_137 = buffer.data(nsf + 137);
    const auto *nsf_138 = buffer.data(nsf + 138);
    const auto *nsf_139 = buffer.data(nsf + 139);
    const auto *nsf_140 = buffer.data(nsf + 140);
    const auto *nsf_145 = buffer.data(nsf + 145);
    const auto *nsf_146 = buffer.data(nsf + 146);
    const auto *nsf_147 = buffer.data(nsf + 147);
    const auto *nsf_149 = buffer.data(nsf + 149);
    const auto *nsf_150 = buffer.data(nsf + 150);
    const auto *nsf_153 = buffer.data(nsf + 153);
    const auto *nsf_156 = buffer.data(nsf + 156);
    const auto *nsf_158 = buffer.data(nsf + 158);
    const auto *nsf_159 = buffer.data(nsf + 159);
    const auto *nsf_165 = buffer.data(nsf + 165);
    const auto *nsf_166 = buffer.data(nsf + 166);
    const auto *nsf_167 = buffer.data(nsf + 167);
    const auto *nsf_168 = buffer.data(nsf + 168);
    const auto *nsf_169 = buffer.data(nsf + 169);

    const auto *nsg1_89 = buffer.data(nsg1 + 89);
    const auto *nsg1_90 = buffer.data(nsg1 + 90);
    const auto *nsg1_93 = buffer.data(nsg1 + 93);
    const auto *nsg1_100 = buffer.data(nsg1 + 100);
    const auto *nsg1_135 = buffer.data(nsg1 + 135);
    const auto *nsg1_138 = buffer.data(nsg1 + 138);
    const auto *nsg1_140 = buffer.data(nsg1 + 140);
    const auto *nsg1_149 = buffer.data(nsg1 + 149);
    const auto *nsg1_150 = buffer.data(nsg1 + 150);
    const auto *nsg1_153 = buffer.data(nsg1 + 153);
    const auto *nsg1_160 = buffer.data(nsg1 + 160);

    const auto *osd0_51 = buffer.data(osd0 + 51);
    const auto *osd0_53 = buffer.data(osd0 + 53);
    const auto *osd0_54 = buffer.data(osd0 + 54);
    const auto *osd0_57 = buffer.data(osd0 + 57);
    const auto *osd0_58 = buffer.data(osd0 + 58);
    const auto *osd0_59 = buffer.data(osd0 + 59);
    const auto *osd0_60 = buffer.data(osd0 + 60);
    const auto *osd0_63 = buffer.data(osd0 + 63);
    const auto *osd0_65 = buffer.data(osd0 + 65);
    const auto *osd0_71 = buffer.data(osd0 + 71);
    const auto *osd0_72 = buffer.data(osd0 + 72);
    const auto *osd0_75 = buffer.data(osd0 + 75);
    const auto *osd0_77 = buffer.data(osd0 + 77);
    const auto *osd0_81 = buffer.data(osd0 + 81);
    const auto *osd0_83 = buffer.data(osd0 + 83);
    const auto *osd0_84 = buffer.data(osd0 + 84);
    const auto *osd0_87 = buffer.data(osd0 + 87);
    const auto *osd0_88 = buffer.data(osd0 + 88);
    const auto *osd0_89 = buffer.data(osd0 + 89);
    const auto *osd0_90 = buffer.data(osd0 + 90);
    const auto *osd0_93 = buffer.data(osd0 + 93);
    const auto *osd0_95 = buffer.data(osd0 + 95);
    const auto *osd0_101 = buffer.data(osd0 + 101);

    const auto *osd1_51 = buffer.data(osd1 + 51);
    const auto *osd1_53 = buffer.data(osd1 + 53);
    const auto *osd1_54 = buffer.data(osd1 + 54);
    const auto *osd1_57 = buffer.data(osd1 + 57);
    const auto *osd1_58 = buffer.data(osd1 + 58);
    const auto *osd1_59 = buffer.data(osd1 + 59);
    const auto *osd1_60 = buffer.data(osd1 + 60);
    const auto *osd1_63 = buffer.data(osd1 + 63);
    const auto *osd1_65 = buffer.data(osd1 + 65);
    const auto *osd1_71 = buffer.data(osd1 + 71);
    const auto *osd1_72 = buffer.data(osd1 + 72);
    const auto *osd1_75 = buffer.data(osd1 + 75);
    const auto *osd1_77 = buffer.data(osd1 + 77);
    const auto *osd1_81 = buffer.data(osd1 + 81);
    const auto *osd1_83 = buffer.data(osd1 + 83);
    const auto *osd1_84 = buffer.data(osd1 + 84);
    const auto *osd1_87 = buffer.data(osd1 + 87);
    const auto *osd1_88 = buffer.data(osd1 + 88);
    const auto *osd1_89 = buffer.data(osd1 + 89);
    const auto *osd1_90 = buffer.data(osd1 + 90);
    const auto *osd1_93 = buffer.data(osd1 + 93);
    const auto *osd1_95 = buffer.data(osd1 + 95);
    const auto *osd1_101 = buffer.data(osd1 + 101);

    const auto *osf_86 = buffer.data(osf + 86);
    const auto *osf_88 = buffer.data(osf + 88);
    const auto *osf_89 = buffer.data(osf + 89);
    const auto *osf_90 = buffer.data(osf + 90);
    const auto *osf_91 = buffer.data(osf + 91);
    const auto *osf_92 = buffer.data(osf + 92);
    const auto *osf_95 = buffer.data(osf + 95);
    const auto *osf_96 = buffer.data(osf + 96);
    const auto *osf_97 = buffer.data(osf + 97);
    const auto *osf_98 = buffer.data(osf + 98);
    const auto *osf_99 = buffer.data(osf + 99);
    const auto *osf_100 = buffer.data(osf + 100);
    const auto *osf_101 = buffer.data(osf + 101);
    const auto *osf_102 = buffer.data(osf + 102);
    const auto *osf_103 = buffer.data(osf + 103);
    const auto *osf_106 = buffer.data(osf + 106);
    const auto *osf_107 = buffer.data(osf + 107);
    const auto *osf_108 = buffer.data(osf + 108);
    const auto *osf_109 = buffer.data(osf + 109);
    const auto *osf_110 = buffer.data(osf + 110);
    const auto *osf_112 = buffer.data(osf + 112);
    const auto *osf_115 = buffer.data(osf + 115);
    const auto *osf_116 = buffer.data(osf + 116);
    const auto *osf_117 = buffer.data(osf + 117);
    const auto *osf_118 = buffer.data(osf + 118);
    const auto *osf_119 = buffer.data(osf + 119);
    const auto *osf_120 = buffer.data(osf + 120);
    const auto *osf_122 = buffer.data(osf + 122);
    const auto *osf_123 = buffer.data(osf + 123);
    const auto *osf_125 = buffer.data(osf + 125);
    const auto *osf_126 = buffer.data(osf + 126);
    const auto *osf_127 = buffer.data(osf + 127);
    const auto *osf_128 = buffer.data(osf + 128);
    const auto *osf_129 = buffer.data(osf + 129);
    const auto *osf_130 = buffer.data(osf + 130);
    const auto *osf_132 = buffer.data(osf + 132);
    const auto *osf_136 = buffer.data(osf + 136);
    const auto *osf_137 = buffer.data(osf + 137);
    const auto *osf_138 = buffer.data(osf + 138);
    const auto *osf_139 = buffer.data(osf + 139);
    const auto *osf_140 = buffer.data(osf + 140);
    const auto *osf_141 = buffer.data(osf + 141);
    const auto *osf_142 = buffer.data(osf + 142);
    const auto *osf_145 = buffer.data(osf + 145);
    const auto *osf_146 = buffer.data(osf + 146);
    const auto *osf_147 = buffer.data(osf + 147);
    const auto *osf_148 = buffer.data(osf + 148);
    const auto *osf_149 = buffer.data(osf + 149);
    const auto *osf_150 = buffer.data(osf + 150);
    const auto *osf_151 = buffer.data(osf + 151);
    const auto *osf_152 = buffer.data(osf + 152);
    const auto *osf_153 = buffer.data(osf + 153);
    const auto *osf_156 = buffer.data(osf + 156);
    const auto *osf_157 = buffer.data(osf + 157);
    const auto *osf_158 = buffer.data(osf + 158);
    const auto *osf_159 = buffer.data(osf + 159);
    const auto *osf_160 = buffer.data(osf + 160);
    const auto *osf_162 = buffer.data(osf + 162);
    const auto *osf_165 = buffer.data(osf + 165);
    const auto *osf_166 = buffer.data(osf + 166);
    const auto *osf_167 = buffer.data(osf + 167);
    const auto *osf_168 = buffer.data(osf + 168);
    const auto *osf_169 = buffer.data(osf + 169);

#pragma omp simd aligned(t_130, t_131, t_132, pc_y, pc_z, nsf_46, nsf_56, nsf_58, osd0_51, \
                         osd0_53, osd1_51, osd1_53, osf_86, osf_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_7 * nsf_56[k]
                   + f_1 * osd0_51[k]
                   - f_2 * osd1_51[k]
                   + f_3 * pc_y[k] * osf_86[k];

        t_131[k] = f_8 * nsf_46[k]
                   + f_3 * pc_z[k] * osf_86[k];

        t_132[k] = f_7 * nsf_58[k]
                   + f_4 * osd0_53[k]
                   - f_5 * osd1_53[k]
                   + f_3 * pc_y[k] * osf_88[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pa_y, pc_x, pc_y, nsg0_89, nsf_59, \
                         nsf_90, nsg1_89, osd0_54, osd1_54, osf_89, \
                         osf_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_7 * nsf_59[k]
                   + f_3 * pc_y[k] * osf_89[k];

        t_134[k] = pa_y[k] * nsg0_89[k]
                   - f_6 * pc_y[k] * nsg1_89[k];

        t_135[k] = f_13 * nsf_90[k]
                   + f_1 * osd0_54[k]
                   - f_2 * osd1_54[k]
                   + f_3 * pc_x[k] * osf_90[k];

        t_136[k] = f_3 * pc_y[k] * osf_90[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pc_y, pc_z, nsf_50, osd0_54, osd1_54, osf_90, \
                         osf_91, osf_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_14 * nsf_50[k]
                   + f_3 * pc_z[k] * osf_90[k];

        t_138[k] = f_4 * osd0_54[k]
                   - f_5 * osd1_54[k]
                   + f_3 * pc_y[k] * osf_91[k];

        t_139[k] = f_3 * pc_y[k] * osf_92[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pc_x, pc_y, nsf_95, nsf_96, nsf_97, \
                         osd0_59, osd1_59, osf_95, osf_96, osf_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_13 * nsf_95[k]
                   + f_4 * osd0_59[k]
                   - f_5 * osd1_59[k]
                   + f_3 * pc_x[k] * osf_95[k];

        t_141[k] = f_13 * nsf_96[k]
                   + f_3 * pc_x[k] * osf_96[k];

        t_142[k] = f_13 * nsf_97[k]
                   + f_3 * pc_x[k] * osf_97[k];

        t_143[k] = f_3 * pc_y[k] * osf_95[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, pc_x, pc_y, nsf_99, osd0_57, osd0_58, osd1_57, \
                         osd1_58, osf_96, osf_97, osf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_13 * nsf_99[k]
                   + f_3 * pc_x[k] * osf_99[k];

        t_145[k] = f_1 * osd0_57[k]
                   - f_2 * osd1_57[k]
                   + f_3 * pc_y[k] * osf_96[k];

        t_146[k] = f_10 * osd0_58[k]
                   - f_11 * osd1_58[k]
                   + f_3 * pc_y[k] * osf_97[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, pc_x, pc_y, pc_z, nsf_59, nsf_100, \
                         osd0_59, osd0_60, osd1_59, osd1_60, osf_98, osf_99, \
                         osf_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_4 * osd0_59[k]
                   - f_5 * osd1_59[k]
                   + f_3 * pc_y[k] * osf_98[k];

        t_148[k] = f_3 * pc_y[k] * osf_99[k];

        t_149[k] = f_14 * nsf_59[k]
                   + f_1 * osd0_59[k]
                   - f_2 * osd1_59[k]
                   + f_3 * pc_z[k] * osf_99[k];

        t_150[k] = f_15 * nsf_100[k]
                   + f_1 * osd0_60[k]
                   - f_2 * osd1_60[k]
                   + f_3 * pc_x[k] * osf_100[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pc_x, pc_y, pc_z, nsf_60, nsf_103, \
                         osd0_63, osd1_63, osf_100, osf_101, osf_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_16 * nsf_60[k]
                   + f_3 * pc_y[k] * osf_100[k];

        t_152[k] = f_3 * pc_z[k] * osf_100[k];

        t_153[k] = f_15 * nsf_103[k]
                   + f_4 * osd0_63[k]
                   - f_5 * osd1_63[k]
                   + f_3 * pc_x[k] * osf_103[k];

        t_154[k] = f_3 * pc_z[k] * osf_101[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pc_x, pc_z, nsf_106, nsf_108, osd0_60, \
                         osd1_60, osf_102, osf_103, osf_106, osf_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_4 * osd0_60[k]
                   - f_5 * osd1_60[k]
                   + f_3 * pc_z[k] * osf_102[k];

        t_156[k] = f_15 * nsf_106[k]
                   + f_3 * pc_x[k] * osf_106[k];

        t_157[k] = f_3 * pc_z[k] * osf_103[k];

        t_158[k] = f_15 * nsf_108[k]
                   + f_3 * pc_x[k] * osf_108[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, t_163, pc_x, pc_y, pc_z, nsf_66, nsf_69, \
                         nsf_109, osd0_63, osd1_63, osf_106, osf_107, \
                         osf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_15 * nsf_109[k]
                   + f_3 * pc_x[k] * osf_109[k];

        t_160[k] = f_16 * nsf_66[k]
                   + f_1 * osd0_63[k]
                   - f_2 * osd1_63[k]
                   + f_3 * pc_y[k] * osf_106[k];

        t_161[k] = f_3 * pc_z[k] * osf_106[k];

        t_162[k] = f_4 * osd0_63[k]
                   - f_5 * osd1_63[k]
                   + f_3 * pc_z[k] * osf_107[k];

        t_163[k] = f_16 * nsf_69[k]
                   + f_3 * pc_y[k] * osf_109[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pa_z, pc_y, pc_z, nsg0_90, nsf_60, \
                         nsf_70, nsg1_90, osd0_65, osd1_65, osf_109, \
                         osf_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_1 * osd0_65[k]
                   - f_2 * osd1_65[k]
                   + f_3 * pc_z[k] * osf_109[k];

        t_165[k] = pa_z[k] * nsg0_90[k]
                   - f_6 * pc_z[k] * nsg1_90[k];

        t_166[k] = f_14 * nsf_70[k]
                   + f_3 * pc_y[k] * osf_110[k];

        t_167[k] = f_7 * nsf_60[k]
                   + f_3 * pc_z[k] * osf_110[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, pa_z, pc_x, pc_y, pc_z, nsg0_93, nsf_72, \
                         nsf_115, nsg1_93, osd0_71, osd1_71, osf_112, \
                         osf_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = pa_z[k] * nsg0_93[k]
                   - f_6 * pc_z[k] * nsg1_93[k];

        t_169[k] = f_14 * nsf_72[k]
                   + f_3 * pc_y[k] * osf_112[k];

        t_170[k] = f_15 * nsf_115[k]
                   + f_4 * osd0_71[k]
                   - f_5 * osd1_71[k]
                   + f_3 * pc_x[k] * osf_115[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, pc_x, nsf_116, nsf_117, nsf_118, nsf_119, \
                         osf_116, osf_117, osf_118, osf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_15 * nsf_116[k]
                   + f_3 * pc_x[k] * osf_116[k];

        t_172[k] = f_15 * nsf_117[k]
                   + f_3 * pc_x[k] * osf_117[k];

        t_173[k] = f_15 * nsf_118[k]
                   + f_3 * pc_x[k] * osf_118[k];

        t_174[k] = f_15 * nsf_119[k]
                   + f_3 * pc_x[k] * osf_119[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, pa_z, pc_y, pc_z, nsg0_100, nsf_66, nsf_78, \
                         nsg1_100, osd0_71, osd1_71, osf_116, osf_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = pa_z[k] * nsg0_100[k]
                   - f_6 * pc_z[k] * nsg1_100[k];

        t_176[k] = f_7 * nsf_66[k]
                   + f_3 * pc_z[k] * osf_116[k];

        t_177[k] = f_14 * nsf_78[k]
                   + f_4 * osd0_71[k]
                   - f_5 * osd1_71[k]
                   + f_3 * pc_y[k] * osf_118[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pc_x, pc_y, pc_z, nsf_69, nsf_79, nsf_120, \
                         osd0_71, osd0_72, osd1_71, osd1_72, osf_119, \
                         osf_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_14 * nsf_79[k]
                   + f_3 * pc_y[k] * osf_119[k];

        t_179[k] = f_7 * nsf_69[k]
                   + f_1 * osd0_71[k]
                   - f_2 * osd1_71[k]
                   + f_3 * pc_z[k] * osf_119[k];

        t_180[k] = f_15 * nsf_120[k]
                   + f_1 * osd0_72[k]
                   - f_2 * osd1_72[k]
                   + f_3 * pc_x[k] * osf_120[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pc_x, pc_y, pc_z, nsf_70, nsf_80, nsf_82, \
                         nsf_123, osd0_75, osd1_75, osf_120, osf_122, \
                         osf_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_8 * nsf_80[k]
                   + f_3 * pc_y[k] * osf_120[k];

        t_182[k] = f_8 * nsf_70[k]
                   + f_3 * pc_z[k] * osf_120[k];

        t_183[k] = f_15 * nsf_123[k]
                   + f_4 * osd0_75[k]
                   - f_5 * osd1_75[k]
                   + f_3 * pc_x[k] * osf_123[k];

        t_184[k] = f_8 * nsf_82[k]
                   + f_3 * pc_y[k] * osf_122[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, nsf_125, nsf_126, nsf_127, nsf_128, \
                         osd0_77, osd1_77, osf_125, osf_126, osf_127, \
                         osf_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_15 * nsf_125[k]
                   + f_4 * osd0_77[k]
                   - f_5 * osd1_77[k]
                   + f_3 * pc_x[k] * osf_125[k];

        t_186[k] = f_15 * nsf_126[k]
                   + f_3 * pc_x[k] * osf_126[k];

        t_187[k] = f_15 * nsf_127[k]
                   + f_3 * pc_x[k] * osf_127[k];

        t_188[k] = f_15 * nsf_128[k]
                   + f_3 * pc_x[k] * osf_128[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pc_x, pc_y, pc_z, nsf_76, nsf_86, nsf_129, \
                         osd0_75, osd1_75, osf_126, osf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_15 * nsf_129[k]
                   + f_3 * pc_x[k] * osf_129[k];

        t_190[k] = f_8 * nsf_86[k]
                   + f_1 * osd0_75[k]
                   - f_2 * osd1_75[k]
                   + f_3 * pc_y[k] * osf_126[k];

        t_191[k] = f_8 * nsf_76[k]
                   + f_3 * pc_z[k] * osf_126[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pa_y, pc_y, pc_z, nsg0_135, nsf_79, \
                         nsf_88, nsf_89, nsg1_135, osd0_77, osd1_77, osf_128, \
                         osf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_8 * nsf_88[k]
                   + f_4 * osd0_77[k]
                   - f_5 * osd1_77[k]
                   + f_3 * pc_y[k] * osf_128[k];

        t_193[k] = f_8 * nsf_89[k]
                   + f_3 * pc_y[k] * osf_129[k];

        t_194[k] = f_8 * nsf_79[k]
                   + f_1 * osd0_77[k]
                   - f_2 * osd1_77[k]
                   + f_3 * pc_z[k] * osf_129[k];

        t_195[k] = pa_y[k] * nsg0_135[k]
                   - f_6 * pc_y[k] * nsg1_135[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pa_y, pc_y, pc_z, nsg0_138, nsf_80, \
                         nsf_90, nsf_91, nsf_92, nsg1_138, osf_130, \
                         osf_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = f_7 * nsf_90[k]
                   + f_3 * pc_y[k] * osf_130[k];

        t_197[k] = f_14 * nsf_80[k]
                   + f_3 * pc_z[k] * osf_130[k];

        t_198[k] = pa_y[k] * nsg0_138[k]
                   + f_8 * nsf_91[k]
                   - f_6 * pc_y[k] * nsg1_138[k];

        t_199[k] = f_7 * nsf_92[k]
                   + f_3 * pc_y[k] * osf_132[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pa_y, pc_x, pc_y, nsg0_140, nsf_136, \
                         nsf_137, nsf_138, nsg1_140, osf_136, osf_137, \
                         osf_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = pa_y[k] * nsg0_140[k]
                   - f_6 * pc_y[k] * nsg1_140[k];

        t_201[k] = f_15 * nsf_136[k]
                   + f_3 * pc_x[k] * osf_136[k];

        t_202[k] = f_15 * nsf_137[k]
                   + f_3 * pc_x[k] * osf_137[k];

        t_203[k] = f_15 * nsf_138[k]
                   + f_3 * pc_x[k] * osf_138[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, pc_x, pc_y, pc_z, nsf_86, nsf_96, nsf_139, \
                         osd0_81, osd1_81, osf_136, osf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_15 * nsf_139[k]
                   + f_3 * pc_x[k] * osf_139[k];

        t_205[k] = f_7 * nsf_96[k]
                   + f_1 * osd0_81[k]
                   - f_2 * osd1_81[k]
                   + f_3 * pc_y[k] * osf_136[k];

        t_206[k] = f_14 * nsf_86[k]
                   + f_3 * pc_z[k] * osf_136[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pa_y, pc_y, nsg0_149, nsf_98, nsf_99, nsg1_149, \
                         osd0_83, osd1_83, osf_138, osf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = f_7 * nsf_98[k]
                   + f_4 * osd0_83[k]
                   - f_5 * osd1_83[k]
                   + f_3 * pc_y[k] * osf_138[k];

        t_208[k] = f_7 * nsf_99[k]
                   + f_3 * pc_y[k] * osf_139[k];

        t_209[k] = pa_y[k] * nsg0_149[k]
                   - f_6 * pc_y[k] * nsg1_149[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, t_214, pc_x, pc_y, pc_z, nsf_90, nsf_140, \
                         osd0_84, osd1_84, osf_140, osf_141, osf_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_15 * nsf_140[k]
                   + f_1 * osd0_84[k]
                   - f_2 * osd1_84[k]
                   + f_3 * pc_x[k] * osf_140[k];

        t_211[k] = f_3 * pc_y[k] * osf_140[k];

        t_212[k] = f_16 * nsf_90[k]
                   + f_3 * pc_z[k] * osf_140[k];

        t_213[k] = f_4 * osd0_84[k]
                   - f_5 * osd1_84[k]
                   + f_3 * pc_y[k] * osf_141[k];

        t_214[k] = f_3 * pc_y[k] * osf_142[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, pc_x, pc_y, nsf_145, nsf_146, nsf_147, \
                         osd0_89, osd1_89, osf_145, osf_146, osf_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_15 * nsf_145[k]
                   + f_4 * osd0_89[k]
                   - f_5 * osd1_89[k]
                   + f_3 * pc_x[k] * osf_145[k];

        t_216[k] = f_15 * nsf_146[k]
                   + f_3 * pc_x[k] * osf_146[k];

        t_217[k] = f_15 * nsf_147[k]
                   + f_3 * pc_x[k] * osf_147[k];

        t_218[k] = f_3 * pc_y[k] * osf_145[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, pc_x, pc_y, nsf_149, osd0_87, osd0_88, osd1_87, \
                         osd1_88, osf_146, osf_147, osf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_15 * nsf_149[k]
                   + f_3 * pc_x[k] * osf_149[k];

        t_220[k] = f_1 * osd0_87[k]
                   - f_2 * osd1_87[k]
                   + f_3 * pc_y[k] * osf_146[k];

        t_221[k] = f_10 * osd0_88[k]
                   - f_11 * osd1_88[k]
                   + f_3 * pc_y[k] * osf_147[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pc_x, pc_y, pc_z, nsf_99, nsf_150, \
                         osd0_89, osd0_90, osd1_89, osd1_90, osf_148, osf_149, \
                         osf_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_4 * osd0_89[k]
                   - f_5 * osd1_89[k]
                   + f_3 * pc_y[k] * osf_148[k];

        t_223[k] = f_3 * pc_y[k] * osf_149[k];

        t_224[k] = f_16 * nsf_99[k]
                   + f_1 * osd0_89[k]
                   - f_2 * osd1_89[k]
                   + f_3 * pc_z[k] * osf_149[k];

        t_225[k] = f_17 * nsf_150[k]
                   + f_1 * osd0_90[k]
                   - f_2 * osd1_90[k]
                   + f_3 * pc_x[k] * osf_150[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, t_229, pc_x, pc_y, pc_z, nsf_100, nsf_153, \
                         osd0_93, osd1_93, osf_150, osf_151, osf_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_18 * nsf_100[k]
                   + f_3 * pc_y[k] * osf_150[k];

        t_227[k] = f_3 * pc_z[k] * osf_150[k];

        t_228[k] = f_17 * nsf_153[k]
                   + f_4 * osd0_93[k]
                   - f_5 * osd1_93[k]
                   + f_3 * pc_x[k] * osf_153[k];

        t_229[k] = f_3 * pc_z[k] * osf_151[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pc_x, pc_z, nsf_156, nsf_158, osd0_90, \
                         osd1_90, osf_152, osf_153, osf_156, osf_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_4 * osd0_90[k]
                   - f_5 * osd1_90[k]
                   + f_3 * pc_z[k] * osf_152[k];

        t_231[k] = f_17 * nsf_156[k]
                   + f_3 * pc_x[k] * osf_156[k];

        t_232[k] = f_3 * pc_z[k] * osf_153[k];

        t_233[k] = f_17 * nsf_158[k]
                   + f_3 * pc_x[k] * osf_158[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, t_238, pc_x, pc_y, pc_z, nsf_106, \
                         nsf_109, nsf_159, osd0_93, osd1_93, osf_156, osf_157, \
                         osf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_17 * nsf_159[k]
                   + f_3 * pc_x[k] * osf_159[k];

        t_235[k] = f_18 * nsf_106[k]
                   + f_1 * osd0_93[k]
                   - f_2 * osd1_93[k]
                   + f_3 * pc_y[k] * osf_156[k];

        t_236[k] = f_3 * pc_z[k] * osf_156[k];

        t_237[k] = f_4 * osd0_93[k]
                   - f_5 * osd1_93[k]
                   + f_3 * pc_z[k] * osf_157[k];

        t_238[k] = f_18 * nsf_109[k]
                   + f_3 * pc_y[k] * osf_159[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, pa_z, pc_y, pc_z, nsg0_150, nsf_100, \
                         nsf_110, nsg1_150, osd0_95, osd1_95, osf_159, \
                         osf_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = f_1 * osd0_95[k]
                   - f_2 * osd1_95[k]
                   + f_3 * pc_z[k] * osf_159[k];

        t_240[k] = pa_z[k] * nsg0_150[k]
                   - f_6 * pc_z[k] * nsg1_150[k];

        t_241[k] = f_16 * nsf_110[k]
                   + f_3 * pc_y[k] * osf_160[k];

        t_242[k] = f_7 * nsf_100[k]
                   + f_3 * pc_z[k] * osf_160[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, pa_z, pc_x, pc_y, pc_z, nsg0_153, nsf_112, \
                         nsf_165, nsg1_153, osd0_101, osd1_101, osf_162, \
                         osf_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = pa_z[k] * nsg0_153[k]
                   - f_6 * pc_z[k] * nsg1_153[k];

        t_244[k] = f_16 * nsf_112[k]
                   + f_3 * pc_y[k] * osf_162[k];

        t_245[k] = f_17 * nsf_165[k]
                   + f_4 * osd0_101[k]
                   - f_5 * osd1_101[k]
                   + f_3 * pc_x[k] * osf_165[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, t_249, pc_x, nsf_166, nsf_167, nsf_168, nsf_169, \
                         osf_166, osf_167, osf_168, osf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = f_17 * nsf_166[k]
                   + f_3 * pc_x[k] * osf_166[k];

        t_247[k] = f_17 * nsf_167[k]
                   + f_3 * pc_x[k] * osf_167[k];

        t_248[k] = f_17 * nsf_168[k]
                   + f_3 * pc_x[k] * osf_168[k];

        t_249[k] = f_17 * nsf_169[k]
                   + f_3 * pc_x[k] * osf_169[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, pa_z, pc_y, pc_z, nsg0_160, nsf_106, nsf_118, \
                         nsg1_160, osd0_101, osd1_101, osf_166, \
                         osf_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = pa_z[k] * nsg0_160[k]
                   - f_6 * pc_z[k] * nsg1_160[k];

        t_251[k] = f_7 * nsf_106[k]
                   + f_3 * pc_z[k] * osf_166[k];

        t_252[k] = f_16 * nsf_118[k]
                   + f_4 * osd0_101[k]
                   - f_5 * osd1_101[k]
                   + f_3 * pc_y[k] * osf_168[k];
    }
}

static auto
compute_prim_osg_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsg0,
                                                          const size_t nsf, const size_t nsg1,
                                                          const size_t osd0, const size_t osd1,
                                                          const size_t osf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_14 = 1.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.5 / q;

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
    auto *t_366 = buffer.data(target + 366);
    auto *t_367 = buffer.data(target + 367);
    auto *t_368 = buffer.data(target + 368);
    auto *t_369 = buffer.data(target + 369);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsg0_210 = buffer.data(nsg0 + 210);
    const auto *nsg0_213 = buffer.data(nsg0 + 213);
    const auto *nsg0_215 = buffer.data(nsg0 + 215);
    const auto *nsg0_224 = buffer.data(nsg0 + 224);
    const auto *nsg0_225 = buffer.data(nsg0 + 225);
    const auto *nsg0_228 = buffer.data(nsg0 + 228);
    const auto *nsg0_235 = buffer.data(nsg0 + 235);

    const auto *nsf_109 = buffer.data(nsf + 109);
    const auto *nsf_110 = buffer.data(nsf + 110);
    const auto *nsf_116 = buffer.data(nsf + 116);
    const auto *nsf_119 = buffer.data(nsf + 119);
    const auto *nsf_120 = buffer.data(nsf + 120);
    const auto *nsf_122 = buffer.data(nsf + 122);
    const auto *nsf_126 = buffer.data(nsf + 126);
    const auto *nsf_128 = buffer.data(nsf + 128);
    const auto *nsf_129 = buffer.data(nsf + 129);
    const auto *nsf_130 = buffer.data(nsf + 130);
    const auto *nsf_132 = buffer.data(nsf + 132);
    const auto *nsf_136 = buffer.data(nsf + 136);
    const auto *nsf_138 = buffer.data(nsf + 138);
    const auto *nsf_139 = buffer.data(nsf + 139);
    const auto *nsf_140 = buffer.data(nsf + 140);
    const auto *nsf_141 = buffer.data(nsf + 141);
    const auto *nsf_142 = buffer.data(nsf + 142);
    const auto *nsf_146 = buffer.data(nsf + 146);
    const auto *nsf_148 = buffer.data(nsf + 148);
    const auto *nsf_149 = buffer.data(nsf + 149);
    const auto *nsf_150 = buffer.data(nsf + 150);
    const auto *nsf_156 = buffer.data(nsf + 156);
    const auto *nsf_159 = buffer.data(nsf + 159);
    const auto *nsf_160 = buffer.data(nsf + 160);
    const auto *nsf_162 = buffer.data(nsf + 162);
    const auto *nsf_166 = buffer.data(nsf + 166);
    const auto *nsf_168 = buffer.data(nsf + 168);
    const auto *nsf_169 = buffer.data(nsf + 169);
    const auto *nsf_170 = buffer.data(nsf + 170);
    const auto *nsf_172 = buffer.data(nsf + 172);
    const auto *nsf_173 = buffer.data(nsf + 173);
    const auto *nsf_175 = buffer.data(nsf + 175);
    const auto *nsf_176 = buffer.data(nsf + 176);
    const auto *nsf_177 = buffer.data(nsf + 177);
    const auto *nsf_178 = buffer.data(nsf + 178);
    const auto *nsf_179 = buffer.data(nsf + 179);
    const auto *nsf_180 = buffer.data(nsf + 180);
    const auto *nsf_182 = buffer.data(nsf + 182);
    const auto *nsf_183 = buffer.data(nsf + 183);
    const auto *nsf_185 = buffer.data(nsf + 185);
    const auto *nsf_186 = buffer.data(nsf + 186);
    const auto *nsf_187 = buffer.data(nsf + 187);
    const auto *nsf_188 = buffer.data(nsf + 188);
    const auto *nsf_189 = buffer.data(nsf + 189);
    const auto *nsf_196 = buffer.data(nsf + 196);
    const auto *nsf_197 = buffer.data(nsf + 197);
    const auto *nsf_198 = buffer.data(nsf + 198);
    const auto *nsf_199 = buffer.data(nsf + 199);
    const auto *nsf_200 = buffer.data(nsf + 200);
    const auto *nsf_205 = buffer.data(nsf + 205);
    const auto *nsf_206 = buffer.data(nsf + 206);
    const auto *nsf_207 = buffer.data(nsf + 207);
    const auto *nsf_209 = buffer.data(nsf + 209);
    const auto *nsf_210 = buffer.data(nsf + 210);
    const auto *nsf_213 = buffer.data(nsf + 213);
    const auto *nsf_216 = buffer.data(nsf + 216);
    const auto *nsf_218 = buffer.data(nsf + 218);
    const auto *nsf_219 = buffer.data(nsf + 219);
    const auto *nsf_225 = buffer.data(nsf + 225);
    const auto *nsf_226 = buffer.data(nsf + 226);
    const auto *nsf_227 = buffer.data(nsf + 227);
    const auto *nsf_228 = buffer.data(nsf + 228);
    const auto *nsf_229 = buffer.data(nsf + 229);
    const auto *nsf_230 = buffer.data(nsf + 230);
    const auto *nsf_233 = buffer.data(nsf + 233);
    const auto *nsf_235 = buffer.data(nsf + 235);
    const auto *nsf_236 = buffer.data(nsf + 236);
    const auto *nsf_237 = buffer.data(nsf + 237);
    const auto *nsf_238 = buffer.data(nsf + 238);
    const auto *nsf_239 = buffer.data(nsf + 239);
    const auto *nsf_240 = buffer.data(nsf + 240);
    const auto *nsf_243 = buffer.data(nsf + 243);
    const auto *nsf_245 = buffer.data(nsf + 245);
    const auto *nsf_246 = buffer.data(nsf + 246);
    const auto *nsf_247 = buffer.data(nsf + 247);
    const auto *nsf_248 = buffer.data(nsf + 248);
    const auto *nsf_249 = buffer.data(nsf + 249);

    const auto *nsg1_210 = buffer.data(nsg1 + 210);
    const auto *nsg1_213 = buffer.data(nsg1 + 213);
    const auto *nsg1_215 = buffer.data(nsg1 + 215);
    const auto *nsg1_224 = buffer.data(nsg1 + 224);
    const auto *nsg1_225 = buffer.data(nsg1 + 225);
    const auto *nsg1_228 = buffer.data(nsg1 + 228);
    const auto *nsg1_235 = buffer.data(nsg1 + 235);

    const auto *osd0_101 = buffer.data(osd0 + 101);
    const auto *osd0_102 = buffer.data(osd0 + 102);
    const auto *osd0_105 = buffer.data(osd0 + 105);
    const auto *osd0_107 = buffer.data(osd0 + 107);
    const auto *osd0_108 = buffer.data(osd0 + 108);
    const auto *osd0_111 = buffer.data(osd0 + 111);
    const auto *osd0_113 = buffer.data(osd0 + 113);
    const auto *osd0_117 = buffer.data(osd0 + 117);
    const auto *osd0_119 = buffer.data(osd0 + 119);
    const auto *osd0_120 = buffer.data(osd0 + 120);
    const auto *osd0_123 = buffer.data(osd0 + 123);
    const auto *osd0_124 = buffer.data(osd0 + 124);
    const auto *osd0_125 = buffer.data(osd0 + 125);
    const auto *osd0_126 = buffer.data(osd0 + 126);
    const auto *osd0_129 = buffer.data(osd0 + 129);
    const auto *osd0_131 = buffer.data(osd0 + 131);
    const auto *osd0_137 = buffer.data(osd0 + 137);
    const auto *osd0_138 = buffer.data(osd0 + 138);
    const auto *osd0_141 = buffer.data(osd0 + 141);
    const auto *osd0_143 = buffer.data(osd0 + 143);
    const auto *osd0_144 = buffer.data(osd0 + 144);
    const auto *osd0_147 = buffer.data(osd0 + 147);
    const auto *osd0_149 = buffer.data(osd0 + 149);

    const auto *osd1_101 = buffer.data(osd1 + 101);
    const auto *osd1_102 = buffer.data(osd1 + 102);
    const auto *osd1_105 = buffer.data(osd1 + 105);
    const auto *osd1_107 = buffer.data(osd1 + 107);
    const auto *osd1_108 = buffer.data(osd1 + 108);
    const auto *osd1_111 = buffer.data(osd1 + 111);
    const auto *osd1_113 = buffer.data(osd1 + 113);
    const auto *osd1_117 = buffer.data(osd1 + 117);
    const auto *osd1_119 = buffer.data(osd1 + 119);
    const auto *osd1_120 = buffer.data(osd1 + 120);
    const auto *osd1_123 = buffer.data(osd1 + 123);
    const auto *osd1_124 = buffer.data(osd1 + 124);
    const auto *osd1_125 = buffer.data(osd1 + 125);
    const auto *osd1_126 = buffer.data(osd1 + 126);
    const auto *osd1_129 = buffer.data(osd1 + 129);
    const auto *osd1_131 = buffer.data(osd1 + 131);
    const auto *osd1_137 = buffer.data(osd1 + 137);
    const auto *osd1_138 = buffer.data(osd1 + 138);
    const auto *osd1_141 = buffer.data(osd1 + 141);
    const auto *osd1_143 = buffer.data(osd1 + 143);
    const auto *osd1_144 = buffer.data(osd1 + 144);
    const auto *osd1_147 = buffer.data(osd1 + 147);
    const auto *osd1_149 = buffer.data(osd1 + 149);

    const auto *osf_169 = buffer.data(osf + 169);
    const auto *osf_170 = buffer.data(osf + 170);
    const auto *osf_172 = buffer.data(osf + 172);
    const auto *osf_173 = buffer.data(osf + 173);
    const auto *osf_175 = buffer.data(osf + 175);
    const auto *osf_176 = buffer.data(osf + 176);
    const auto *osf_177 = buffer.data(osf + 177);
    const auto *osf_178 = buffer.data(osf + 178);
    const auto *osf_179 = buffer.data(osf + 179);
    const auto *osf_180 = buffer.data(osf + 180);
    const auto *osf_182 = buffer.data(osf + 182);
    const auto *osf_183 = buffer.data(osf + 183);
    const auto *osf_185 = buffer.data(osf + 185);
    const auto *osf_186 = buffer.data(osf + 186);
    const auto *osf_187 = buffer.data(osf + 187);
    const auto *osf_188 = buffer.data(osf + 188);
    const auto *osf_189 = buffer.data(osf + 189);
    const auto *osf_190 = buffer.data(osf + 190);
    const auto *osf_192 = buffer.data(osf + 192);
    const auto *osf_196 = buffer.data(osf + 196);
    const auto *osf_197 = buffer.data(osf + 197);
    const auto *osf_198 = buffer.data(osf + 198);
    const auto *osf_199 = buffer.data(osf + 199);
    const auto *osf_200 = buffer.data(osf + 200);
    const auto *osf_201 = buffer.data(osf + 201);
    const auto *osf_202 = buffer.data(osf + 202);
    const auto *osf_205 = buffer.data(osf + 205);
    const auto *osf_206 = buffer.data(osf + 206);
    const auto *osf_207 = buffer.data(osf + 207);
    const auto *osf_208 = buffer.data(osf + 208);
    const auto *osf_209 = buffer.data(osf + 209);
    const auto *osf_210 = buffer.data(osf + 210);
    const auto *osf_211 = buffer.data(osf + 211);
    const auto *osf_212 = buffer.data(osf + 212);
    const auto *osf_213 = buffer.data(osf + 213);
    const auto *osf_216 = buffer.data(osf + 216);
    const auto *osf_217 = buffer.data(osf + 217);
    const auto *osf_218 = buffer.data(osf + 218);
    const auto *osf_219 = buffer.data(osf + 219);
    const auto *osf_220 = buffer.data(osf + 220);
    const auto *osf_222 = buffer.data(osf + 222);
    const auto *osf_225 = buffer.data(osf + 225);
    const auto *osf_226 = buffer.data(osf + 226);
    const auto *osf_227 = buffer.data(osf + 227);
    const auto *osf_228 = buffer.data(osf + 228);
    const auto *osf_229 = buffer.data(osf + 229);
    const auto *osf_230 = buffer.data(osf + 230);
    const auto *osf_232 = buffer.data(osf + 232);
    const auto *osf_233 = buffer.data(osf + 233);
    const auto *osf_235 = buffer.data(osf + 235);
    const auto *osf_236 = buffer.data(osf + 236);
    const auto *osf_237 = buffer.data(osf + 237);
    const auto *osf_238 = buffer.data(osf + 238);
    const auto *osf_239 = buffer.data(osf + 239);
    const auto *osf_240 = buffer.data(osf + 240);
    const auto *osf_242 = buffer.data(osf + 242);
    const auto *osf_243 = buffer.data(osf + 243);
    const auto *osf_245 = buffer.data(osf + 245);
    const auto *osf_246 = buffer.data(osf + 246);
    const auto *osf_247 = buffer.data(osf + 247);
    const auto *osf_248 = buffer.data(osf + 248);
    const auto *osf_249 = buffer.data(osf + 249);

#pragma omp simd aligned(t_253, t_254, t_255, pc_x, pc_y, pc_z, nsf_109, nsf_119, nsf_170, \
                         osd0_101, osd0_102, osd1_101, osd1_102, osf_169, \
                         osf_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_16 * nsf_119[k]
                   + f_3 * pc_y[k] * osf_169[k];

        t_254[k] = f_7 * nsf_109[k]
                   + f_1 * osd0_101[k]
                   - f_2 * osd1_101[k]
                   + f_3 * pc_z[k] * osf_169[k];

        t_255[k] = f_17 * nsf_170[k]
                   + f_1 * osd0_102[k]
                   - f_2 * osd1_102[k]
                   + f_3 * pc_x[k] * osf_170[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, t_259, pc_x, pc_y, pc_z, nsf_110, nsf_120, \
                         nsf_122, nsf_173, osd0_105, osd1_105, osf_170, osf_172, \
                         osf_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_14 * nsf_120[k]
                   + f_3 * pc_y[k] * osf_170[k];

        t_257[k] = f_8 * nsf_110[k]
                   + f_3 * pc_z[k] * osf_170[k];

        t_258[k] = f_17 * nsf_173[k]
                   + f_4 * osd0_105[k]
                   - f_5 * osd1_105[k]
                   + f_3 * pc_x[k] * osf_173[k];

        t_259[k] = f_14 * nsf_122[k]
                   + f_3 * pc_y[k] * osf_172[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pc_x, nsf_175, nsf_176, nsf_177, nsf_178, \
                         osd0_107, osd1_107, osf_175, osf_176, osf_177, \
                         osf_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = f_17 * nsf_175[k]
                   + f_4 * osd0_107[k]
                   - f_5 * osd1_107[k]
                   + f_3 * pc_x[k] * osf_175[k];

        t_261[k] = f_17 * nsf_176[k]
                   + f_3 * pc_x[k] * osf_176[k];

        t_262[k] = f_17 * nsf_177[k]
                   + f_3 * pc_x[k] * osf_177[k];

        t_263[k] = f_17 * nsf_178[k]
                   + f_3 * pc_x[k] * osf_178[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pc_x, pc_y, pc_z, nsf_116, nsf_126, nsf_179, \
                         osd0_105, osd1_105, osf_176, osf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_17 * nsf_179[k]
                   + f_3 * pc_x[k] * osf_179[k];

        t_265[k] = f_14 * nsf_126[k]
                   + f_1 * osd0_105[k]
                   - f_2 * osd1_105[k]
                   + f_3 * pc_y[k] * osf_176[k];

        t_266[k] = f_8 * nsf_116[k]
                   + f_3 * pc_z[k] * osf_176[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, pc_y, pc_z, nsf_119, nsf_128, nsf_129, osd0_107, \
                         osd1_107, osf_178, osf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_14 * nsf_128[k]
                   + f_4 * osd0_107[k]
                   - f_5 * osd1_107[k]
                   + f_3 * pc_y[k] * osf_178[k];

        t_268[k] = f_14 * nsf_129[k]
                   + f_3 * pc_y[k] * osf_179[k];

        t_269[k] = f_8 * nsf_119[k]
                   + f_1 * osd0_107[k]
                   - f_2 * osd1_107[k]
                   + f_3 * pc_z[k] * osf_179[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, pc_x, pc_y, pc_z, nsf_120, nsf_130, nsf_180, \
                         osd0_108, osd1_108, osf_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_17 * nsf_180[k]
                   + f_1 * osd0_108[k]
                   - f_2 * osd1_108[k]
                   + f_3 * pc_x[k] * osf_180[k];

        t_271[k] = f_8 * nsf_130[k]
                   + f_3 * pc_y[k] * osf_180[k];

        t_272[k] = f_14 * nsf_120[k]
                   + f_3 * pc_z[k] * osf_180[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, pc_x, pc_y, nsf_132, nsf_183, nsf_185, osd0_111, \
                         osd0_113, osd1_111, osd1_113, osf_182, osf_183, \
                         osf_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_17 * nsf_183[k]
                   + f_4 * osd0_111[k]
                   - f_5 * osd1_111[k]
                   + f_3 * pc_x[k] * osf_183[k];

        t_274[k] = f_8 * nsf_132[k]
                   + f_3 * pc_y[k] * osf_182[k];

        t_275[k] = f_17 * nsf_185[k]
                   + f_4 * osd0_113[k]
                   - f_5 * osd1_113[k]
                   + f_3 * pc_x[k] * osf_185[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pc_x, nsf_186, nsf_187, nsf_188, nsf_189, \
                         osf_186, osf_187, osf_188, osf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_17 * nsf_186[k]
                   + f_3 * pc_x[k] * osf_186[k];

        t_277[k] = f_17 * nsf_187[k]
                   + f_3 * pc_x[k] * osf_187[k];

        t_278[k] = f_17 * nsf_188[k]
                   + f_3 * pc_x[k] * osf_188[k];

        t_279[k] = f_17 * nsf_189[k]
                   + f_3 * pc_x[k] * osf_189[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, pc_y, pc_z, nsf_126, nsf_136, nsf_138, osd0_111, \
                         osd0_113, osd1_111, osd1_113, osf_186, \
                         osf_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_8 * nsf_136[k]
                   + f_1 * osd0_111[k]
                   - f_2 * osd1_111[k]
                   + f_3 * pc_y[k] * osf_186[k];

        t_281[k] = f_14 * nsf_126[k]
                   + f_3 * pc_z[k] * osf_186[k];

        t_282[k] = f_8 * nsf_138[k]
                   + f_4 * osd0_113[k]
                   - f_5 * osd1_113[k]
                   + f_3 * pc_y[k] * osf_188[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, pa_y, pc_y, pc_z, nsg0_210, nsf_129, \
                         nsf_139, nsf_140, nsg1_210, osd0_113, osd1_113, osf_189, \
                         osf_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_8 * nsf_139[k]
                   + f_3 * pc_y[k] * osf_189[k];

        t_284[k] = f_14 * nsf_129[k]
                   + f_1 * osd0_113[k]
                   - f_2 * osd1_113[k]
                   + f_3 * pc_z[k] * osf_189[k];

        t_285[k] = pa_y[k] * nsg0_210[k]
                   - f_6 * pc_y[k] * nsg1_210[k];

        t_286[k] = f_7 * nsf_140[k]
                   + f_3 * pc_y[k] * osf_190[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, t_290, pa_y, pc_y, pc_z, nsg0_213, nsg0_215, \
                         nsf_130, nsf_141, nsf_142, nsg1_213, nsg1_215, osf_190, \
                         osf_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = f_16 * nsf_130[k]
                   + f_3 * pc_z[k] * osf_190[k];

        t_288[k] = pa_y[k] * nsg0_213[k]
                   + f_8 * nsf_141[k]
                   - f_6 * pc_y[k] * nsg1_213[k];

        t_289[k] = f_7 * nsf_142[k]
                   + f_3 * pc_y[k] * osf_192[k];

        t_290[k] = pa_y[k] * nsg0_215[k]
                   - f_6 * pc_y[k] * nsg1_215[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, t_294, pc_x, nsf_196, nsf_197, nsf_198, nsf_199, \
                         osf_196, osf_197, osf_198, osf_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_17 * nsf_196[k]
                   + f_3 * pc_x[k] * osf_196[k];

        t_292[k] = f_17 * nsf_197[k]
                   + f_3 * pc_x[k] * osf_197[k];

        t_293[k] = f_17 * nsf_198[k]
                   + f_3 * pc_x[k] * osf_198[k];

        t_294[k] = f_17 * nsf_199[k]
                   + f_3 * pc_x[k] * osf_199[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, pc_y, pc_z, nsf_136, nsf_146, nsf_148, osd0_117, \
                         osd0_119, osd1_117, osd1_119, osf_196, \
                         osf_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = f_7 * nsf_146[k]
                   + f_1 * osd0_117[k]
                   - f_2 * osd1_117[k]
                   + f_3 * pc_y[k] * osf_196[k];

        t_296[k] = f_16 * nsf_136[k]
                   + f_3 * pc_z[k] * osf_196[k];

        t_297[k] = f_7 * nsf_148[k]
                   + f_4 * osd0_119[k]
                   - f_5 * osd1_119[k]
                   + f_3 * pc_y[k] * osf_198[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, pa_y, pc_x, pc_y, nsg0_224, nsf_149, \
                         nsf_200, nsg1_224, osd0_120, osd1_120, osf_199, \
                         osf_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_7 * nsf_149[k]
                   + f_3 * pc_y[k] * osf_199[k];

        t_299[k] = pa_y[k] * nsg0_224[k]
                   - f_6 * pc_y[k] * nsg1_224[k];

        t_300[k] = f_17 * nsf_200[k]
                   + f_1 * osd0_120[k]
                   - f_2 * osd1_120[k]
                   + f_3 * pc_x[k] * osf_200[k];

        t_301[k] = f_3 * pc_y[k] * osf_200[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, pc_y, pc_z, nsf_140, osd0_120, osd1_120, \
                         osf_200, osf_201, osf_202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = f_18 * nsf_140[k]
                   + f_3 * pc_z[k] * osf_200[k];

        t_303[k] = f_4 * osd0_120[k]
                   - f_5 * osd1_120[k]
                   + f_3 * pc_y[k] * osf_201[k];

        t_304[k] = f_3 * pc_y[k] * osf_202[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pc_x, pc_y, nsf_205, nsf_206, nsf_207, \
                         osd0_125, osd1_125, osf_205, osf_206, \
                         osf_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_17 * nsf_205[k]
                   + f_4 * osd0_125[k]
                   - f_5 * osd1_125[k]
                   + f_3 * pc_x[k] * osf_205[k];

        t_306[k] = f_17 * nsf_206[k]
                   + f_3 * pc_x[k] * osf_206[k];

        t_307[k] = f_17 * nsf_207[k]
                   + f_3 * pc_x[k] * osf_207[k];

        t_308[k] = f_3 * pc_y[k] * osf_205[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, pc_x, pc_y, nsf_209, osd0_123, osd0_124, \
                         osd1_123, osd1_124, osf_206, osf_207, \
                         osf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_17 * nsf_209[k]
                   + f_3 * pc_x[k] * osf_209[k];

        t_310[k] = f_1 * osd0_123[k]
                   - f_2 * osd1_123[k]
                   + f_3 * pc_y[k] * osf_206[k];

        t_311[k] = f_10 * osd0_124[k]
                   - f_11 * osd1_124[k]
                   + f_3 * pc_y[k] * osf_207[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, t_315, pc_x, pc_y, pc_z, nsf_149, nsf_210, \
                         osd0_125, osd0_126, osd1_125, osd1_126, osf_208, osf_209, \
                         osf_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_4 * osd0_125[k]
                   - f_5 * osd1_125[k]
                   + f_3 * pc_y[k] * osf_208[k];

        t_313[k] = f_3 * pc_y[k] * osf_209[k];

        t_314[k] = f_18 * nsf_149[k]
                   + f_1 * osd0_125[k]
                   - f_2 * osd1_125[k]
                   + f_3 * pc_z[k] * osf_209[k];

        t_315[k] = f_18 * nsf_210[k]
                   + f_1 * osd0_126[k]
                   - f_2 * osd1_126[k]
                   + f_3 * pc_x[k] * osf_210[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, t_319, pc_x, pc_y, pc_z, nsf_150, nsf_213, \
                         osd0_129, osd1_129, osf_210, osf_211, \
                         osf_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = f_17 * nsf_150[k]
                   + f_3 * pc_y[k] * osf_210[k];

        t_317[k] = f_3 * pc_z[k] * osf_210[k];

        t_318[k] = f_18 * nsf_213[k]
                   + f_4 * osd0_129[k]
                   - f_5 * osd1_129[k]
                   + f_3 * pc_x[k] * osf_213[k];

        t_319[k] = f_3 * pc_z[k] * osf_211[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, pc_x, pc_z, nsf_216, nsf_218, osd0_126, \
                         osd1_126, osf_212, osf_213, osf_216, osf_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_4 * osd0_126[k]
                   - f_5 * osd1_126[k]
                   + f_3 * pc_z[k] * osf_212[k];

        t_321[k] = f_18 * nsf_216[k]
                   + f_3 * pc_x[k] * osf_216[k];

        t_322[k] = f_3 * pc_z[k] * osf_213[k];

        t_323[k] = f_18 * nsf_218[k]
                   + f_3 * pc_x[k] * osf_218[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, t_328, pc_x, pc_y, pc_z, nsf_156, \
                         nsf_159, nsf_219, osd0_129, osd1_129, osf_216, osf_217, \
                         osf_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_18 * nsf_219[k]
                   + f_3 * pc_x[k] * osf_219[k];

        t_325[k] = f_17 * nsf_156[k]
                   + f_1 * osd0_129[k]
                   - f_2 * osd1_129[k]
                   + f_3 * pc_y[k] * osf_216[k];

        t_326[k] = f_3 * pc_z[k] * osf_216[k];

        t_327[k] = f_4 * osd0_129[k]
                   - f_5 * osd1_129[k]
                   + f_3 * pc_z[k] * osf_217[k];

        t_328[k] = f_17 * nsf_159[k]
                   + f_3 * pc_y[k] * osf_219[k];
    }

#pragma omp simd aligned(t_329, t_330, t_331, t_332, pa_z, pc_y, pc_z, nsg0_225, nsf_150, \
                         nsf_160, nsg1_225, osd0_131, osd1_131, osf_219, \
                         osf_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_329[k] = f_1 * osd0_131[k]
                   - f_2 * osd1_131[k]
                   + f_3 * pc_z[k] * osf_219[k];

        t_330[k] = pa_z[k] * nsg0_225[k]
                   - f_6 * pc_z[k] * nsg1_225[k];

        t_331[k] = f_18 * nsf_160[k]
                   + f_3 * pc_y[k] * osf_220[k];

        t_332[k] = f_7 * nsf_150[k]
                   + f_3 * pc_z[k] * osf_220[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, pa_z, pc_x, pc_y, pc_z, nsg0_228, nsf_162, \
                         nsf_225, nsg1_228, osd0_137, osd1_137, osf_222, \
                         osf_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = pa_z[k] * nsg0_228[k]
                   - f_6 * pc_z[k] * nsg1_228[k];

        t_334[k] = f_18 * nsf_162[k]
                   + f_3 * pc_y[k] * osf_222[k];

        t_335[k] = f_18 * nsf_225[k]
                   + f_4 * osd0_137[k]
                   - f_5 * osd1_137[k]
                   + f_3 * pc_x[k] * osf_225[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, t_339, pc_x, nsf_226, nsf_227, nsf_228, nsf_229, \
                         osf_226, osf_227, osf_228, osf_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = f_18 * nsf_226[k]
                   + f_3 * pc_x[k] * osf_226[k];

        t_337[k] = f_18 * nsf_227[k]
                   + f_3 * pc_x[k] * osf_227[k];

        t_338[k] = f_18 * nsf_228[k]
                   + f_3 * pc_x[k] * osf_228[k];

        t_339[k] = f_18 * nsf_229[k]
                   + f_3 * pc_x[k] * osf_229[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, pa_z, pc_y, pc_z, nsg0_235, nsf_156, nsf_168, \
                         nsg1_235, osd0_137, osd1_137, osf_226, \
                         osf_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = pa_z[k] * nsg0_235[k]
                   - f_6 * pc_z[k] * nsg1_235[k];

        t_341[k] = f_7 * nsf_156[k]
                   + f_3 * pc_z[k] * osf_226[k];

        t_342[k] = f_18 * nsf_168[k]
                   + f_4 * osd0_137[k]
                   - f_5 * osd1_137[k]
                   + f_3 * pc_y[k] * osf_228[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, pc_x, pc_y, pc_z, nsf_159, nsf_169, nsf_230, \
                         osd0_137, osd0_138, osd1_137, osd1_138, osf_229, \
                         osf_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = f_18 * nsf_169[k]
                   + f_3 * pc_y[k] * osf_229[k];

        t_344[k] = f_7 * nsf_159[k]
                   + f_1 * osd0_137[k]
                   - f_2 * osd1_137[k]
                   + f_3 * pc_z[k] * osf_229[k];

        t_345[k] = f_18 * nsf_230[k]
                   + f_1 * osd0_138[k]
                   - f_2 * osd1_138[k]
                   + f_3 * pc_x[k] * osf_230[k];
    }

#pragma omp simd aligned(t_346, t_347, t_348, t_349, pc_x, pc_y, pc_z, nsf_160, nsf_170, \
                         nsf_172, nsf_233, osd0_141, osd1_141, osf_230, osf_232, \
                         osf_233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_346[k] = f_16 * nsf_170[k]
                   + f_3 * pc_y[k] * osf_230[k];

        t_347[k] = f_8 * nsf_160[k]
                   + f_3 * pc_z[k] * osf_230[k];

        t_348[k] = f_18 * nsf_233[k]
                   + f_4 * osd0_141[k]
                   - f_5 * osd1_141[k]
                   + f_3 * pc_x[k] * osf_233[k];

        t_349[k] = f_16 * nsf_172[k]
                   + f_3 * pc_y[k] * osf_232[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, pc_x, nsf_235, nsf_236, nsf_237, nsf_238, \
                         osd0_143, osd1_143, osf_235, osf_236, osf_237, \
                         osf_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = f_18 * nsf_235[k]
                   + f_4 * osd0_143[k]
                   - f_5 * osd1_143[k]
                   + f_3 * pc_x[k] * osf_235[k];

        t_351[k] = f_18 * nsf_236[k]
                   + f_3 * pc_x[k] * osf_236[k];

        t_352[k] = f_18 * nsf_237[k]
                   + f_3 * pc_x[k] * osf_237[k];

        t_353[k] = f_18 * nsf_238[k]
                   + f_3 * pc_x[k] * osf_238[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, pc_x, pc_y, pc_z, nsf_166, nsf_176, nsf_239, \
                         osd0_141, osd1_141, osf_236, osf_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_18 * nsf_239[k]
                   + f_3 * pc_x[k] * osf_239[k];

        t_355[k] = f_16 * nsf_176[k]
                   + f_1 * osd0_141[k]
                   - f_2 * osd1_141[k]
                   + f_3 * pc_y[k] * osf_236[k];

        t_356[k] = f_8 * nsf_166[k]
                   + f_3 * pc_z[k] * osf_236[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, pc_y, pc_z, nsf_169, nsf_178, nsf_179, osd0_143, \
                         osd1_143, osf_238, osf_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = f_16 * nsf_178[k]
                   + f_4 * osd0_143[k]
                   - f_5 * osd1_143[k]
                   + f_3 * pc_y[k] * osf_238[k];

        t_358[k] = f_16 * nsf_179[k]
                   + f_3 * pc_y[k] * osf_239[k];

        t_359[k] = f_8 * nsf_169[k]
                   + f_1 * osd0_143[k]
                   - f_2 * osd1_143[k]
                   + f_3 * pc_z[k] * osf_239[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, pc_x, pc_y, pc_z, nsf_170, nsf_180, nsf_240, \
                         osd0_144, osd1_144, osf_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = f_18 * nsf_240[k]
                   + f_1 * osd0_144[k]
                   - f_2 * osd1_144[k]
                   + f_3 * pc_x[k] * osf_240[k];

        t_361[k] = f_14 * nsf_180[k]
                   + f_3 * pc_y[k] * osf_240[k];

        t_362[k] = f_14 * nsf_170[k]
                   + f_3 * pc_z[k] * osf_240[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, pc_x, pc_y, nsf_182, nsf_243, nsf_245, osd0_147, \
                         osd0_149, osd1_147, osd1_149, osf_242, osf_243, \
                         osf_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = f_18 * nsf_243[k]
                   + f_4 * osd0_147[k]
                   - f_5 * osd1_147[k]
                   + f_3 * pc_x[k] * osf_243[k];

        t_364[k] = f_14 * nsf_182[k]
                   + f_3 * pc_y[k] * osf_242[k];

        t_365[k] = f_18 * nsf_245[k]
                   + f_4 * osd0_149[k]
                   - f_5 * osd1_149[k]
                   + f_3 * pc_x[k] * osf_245[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pc_x, nsf_246, nsf_247, nsf_248, nsf_249, \
                         osf_246, osf_247, osf_248, osf_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_18 * nsf_246[k]
                   + f_3 * pc_x[k] * osf_246[k];

        t_367[k] = f_18 * nsf_247[k]
                   + f_3 * pc_x[k] * osf_247[k];

        t_368[k] = f_18 * nsf_248[k]
                   + f_3 * pc_x[k] * osf_248[k];

        t_369[k] = f_18 * nsf_249[k]
                   + f_3 * pc_x[k] * osf_249[k];
    }
}

static auto
compute_prim_osg_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsg0,
                                                          const size_t nsf, const size_t nsg1,
                                                          const size_t osd0, const size_t osd1,
                                                          const size_t osf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_14 = 1.5 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsg0_300 = buffer.data(nsg0 + 300);
    const auto *nsg0_303 = buffer.data(nsg0 + 303);
    const auto *nsg0_305 = buffer.data(nsg0 + 305);
    const auto *nsg0_314 = buffer.data(nsg0 + 314);
    const auto *nsg0_315 = buffer.data(nsg0 + 315);
    const auto *nsg0_318 = buffer.data(nsg0 + 318);
    const auto *nsg0_325 = buffer.data(nsg0 + 325);

    const auto *nsf_176 = buffer.data(nsf + 176);
    const auto *nsf_179 = buffer.data(nsf + 179);
    const auto *nsf_180 = buffer.data(nsf + 180);
    const auto *nsf_186 = buffer.data(nsf + 186);
    const auto *nsf_188 = buffer.data(nsf + 188);
    const auto *nsf_189 = buffer.data(nsf + 189);
    const auto *nsf_190 = buffer.data(nsf + 190);
    const auto *nsf_192 = buffer.data(nsf + 192);
    const auto *nsf_196 = buffer.data(nsf + 196);
    const auto *nsf_198 = buffer.data(nsf + 198);
    const auto *nsf_199 = buffer.data(nsf + 199);
    const auto *nsf_200 = buffer.data(nsf + 200);
    const auto *nsf_201 = buffer.data(nsf + 201);
    const auto *nsf_202 = buffer.data(nsf + 202);
    const auto *nsf_206 = buffer.data(nsf + 206);
    const auto *nsf_208 = buffer.data(nsf + 208);
    const auto *nsf_209 = buffer.data(nsf + 209);
    const auto *nsf_210 = buffer.data(nsf + 210);
    const auto *nsf_216 = buffer.data(nsf + 216);
    const auto *nsf_219 = buffer.data(nsf + 219);
    const auto *nsf_220 = buffer.data(nsf + 220);
    const auto *nsf_222 = buffer.data(nsf + 222);
    const auto *nsf_226 = buffer.data(nsf + 226);
    const auto *nsf_228 = buffer.data(nsf + 228);
    const auto *nsf_229 = buffer.data(nsf + 229);
    const auto *nsf_230 = buffer.data(nsf + 230);
    const auto *nsf_232 = buffer.data(nsf + 232);
    const auto *nsf_236 = buffer.data(nsf + 236);
    const auto *nsf_238 = buffer.data(nsf + 238);
    const auto *nsf_239 = buffer.data(nsf + 239);
    const auto *nsf_240 = buffer.data(nsf + 240);
    const auto *nsf_242 = buffer.data(nsf + 242);
    const auto *nsf_246 = buffer.data(nsf + 246);
    const auto *nsf_248 = buffer.data(nsf + 248);
    const auto *nsf_249 = buffer.data(nsf + 249);
    const auto *nsf_250 = buffer.data(nsf + 250);
    const auto *nsf_252 = buffer.data(nsf + 252);
    const auto *nsf_253 = buffer.data(nsf + 253);
    const auto *nsf_255 = buffer.data(nsf + 255);
    const auto *nsf_256 = buffer.data(nsf + 256);
    const auto *nsf_257 = buffer.data(nsf + 257);
    const auto *nsf_258 = buffer.data(nsf + 258);
    const auto *nsf_259 = buffer.data(nsf + 259);
    const auto *nsf_266 = buffer.data(nsf + 266);
    const auto *nsf_267 = buffer.data(nsf + 267);
    const auto *nsf_268 = buffer.data(nsf + 268);
    const auto *nsf_269 = buffer.data(nsf + 269);
    const auto *nsf_270 = buffer.data(nsf + 270);
    const auto *nsf_275 = buffer.data(nsf + 275);
    const auto *nsf_276 = buffer.data(nsf + 276);
    const auto *nsf_277 = buffer.data(nsf + 277);
    const auto *nsf_279 = buffer.data(nsf + 279);
    const auto *nsf_280 = buffer.data(nsf + 280);
    const auto *nsf_283 = buffer.data(nsf + 283);
    const auto *nsf_286 = buffer.data(nsf + 286);
    const auto *nsf_288 = buffer.data(nsf + 288);
    const auto *nsf_289 = buffer.data(nsf + 289);
    const auto *nsf_295 = buffer.data(nsf + 295);
    const auto *nsf_296 = buffer.data(nsf + 296);
    const auto *nsf_297 = buffer.data(nsf + 297);
    const auto *nsf_298 = buffer.data(nsf + 298);
    const auto *nsf_299 = buffer.data(nsf + 299);
    const auto *nsf_300 = buffer.data(nsf + 300);
    const auto *nsf_303 = buffer.data(nsf + 303);
    const auto *nsf_305 = buffer.data(nsf + 305);
    const auto *nsf_306 = buffer.data(nsf + 306);
    const auto *nsf_307 = buffer.data(nsf + 307);
    const auto *nsf_308 = buffer.data(nsf + 308);
    const auto *nsf_309 = buffer.data(nsf + 309);
    const auto *nsf_310 = buffer.data(nsf + 310);
    const auto *nsf_313 = buffer.data(nsf + 313);
    const auto *nsf_315 = buffer.data(nsf + 315);
    const auto *nsf_316 = buffer.data(nsf + 316);
    const auto *nsf_317 = buffer.data(nsf + 317);
    const auto *nsf_318 = buffer.data(nsf + 318);
    const auto *nsf_319 = buffer.data(nsf + 319);
    const auto *nsf_320 = buffer.data(nsf + 320);
    const auto *nsf_323 = buffer.data(nsf + 323);

    const auto *nsg1_300 = buffer.data(nsg1 + 300);
    const auto *nsg1_303 = buffer.data(nsg1 + 303);
    const auto *nsg1_305 = buffer.data(nsg1 + 305);
    const auto *nsg1_314 = buffer.data(nsg1 + 314);
    const auto *nsg1_315 = buffer.data(nsg1 + 315);
    const auto *nsg1_318 = buffer.data(nsg1 + 318);
    const auto *nsg1_325 = buffer.data(nsg1 + 325);

    const auto *osd0_147 = buffer.data(osd0 + 147);
    const auto *osd0_149 = buffer.data(osd0 + 149);
    const auto *osd0_150 = buffer.data(osd0 + 150);
    const auto *osd0_153 = buffer.data(osd0 + 153);
    const auto *osd0_155 = buffer.data(osd0 + 155);
    const auto *osd0_159 = buffer.data(osd0 + 159);
    const auto *osd0_161 = buffer.data(osd0 + 161);
    const auto *osd0_162 = buffer.data(osd0 + 162);
    const auto *osd0_165 = buffer.data(osd0 + 165);
    const auto *osd0_166 = buffer.data(osd0 + 166);
    const auto *osd0_167 = buffer.data(osd0 + 167);
    const auto *osd0_168 = buffer.data(osd0 + 168);
    const auto *osd0_171 = buffer.data(osd0 + 171);
    const auto *osd0_173 = buffer.data(osd0 + 173);
    const auto *osd0_179 = buffer.data(osd0 + 179);
    const auto *osd0_180 = buffer.data(osd0 + 180);
    const auto *osd0_183 = buffer.data(osd0 + 183);
    const auto *osd0_185 = buffer.data(osd0 + 185);
    const auto *osd0_186 = buffer.data(osd0 + 186);
    const auto *osd0_189 = buffer.data(osd0 + 189);
    const auto *osd0_191 = buffer.data(osd0 + 191);
    const auto *osd0_192 = buffer.data(osd0 + 192);
    const auto *osd0_195 = buffer.data(osd0 + 195);

    const auto *osd1_147 = buffer.data(osd1 + 147);
    const auto *osd1_149 = buffer.data(osd1 + 149);
    const auto *osd1_150 = buffer.data(osd1 + 150);
    const auto *osd1_153 = buffer.data(osd1 + 153);
    const auto *osd1_155 = buffer.data(osd1 + 155);
    const auto *osd1_159 = buffer.data(osd1 + 159);
    const auto *osd1_161 = buffer.data(osd1 + 161);
    const auto *osd1_162 = buffer.data(osd1 + 162);
    const auto *osd1_165 = buffer.data(osd1 + 165);
    const auto *osd1_166 = buffer.data(osd1 + 166);
    const auto *osd1_167 = buffer.data(osd1 + 167);
    const auto *osd1_168 = buffer.data(osd1 + 168);
    const auto *osd1_171 = buffer.data(osd1 + 171);
    const auto *osd1_173 = buffer.data(osd1 + 173);
    const auto *osd1_179 = buffer.data(osd1 + 179);
    const auto *osd1_180 = buffer.data(osd1 + 180);
    const auto *osd1_183 = buffer.data(osd1 + 183);
    const auto *osd1_185 = buffer.data(osd1 + 185);
    const auto *osd1_186 = buffer.data(osd1 + 186);
    const auto *osd1_189 = buffer.data(osd1 + 189);
    const auto *osd1_191 = buffer.data(osd1 + 191);
    const auto *osd1_192 = buffer.data(osd1 + 192);
    const auto *osd1_195 = buffer.data(osd1 + 195);

    const auto *osf_246 = buffer.data(osf + 246);
    const auto *osf_248 = buffer.data(osf + 248);
    const auto *osf_249 = buffer.data(osf + 249);
    const auto *osf_250 = buffer.data(osf + 250);
    const auto *osf_252 = buffer.data(osf + 252);
    const auto *osf_253 = buffer.data(osf + 253);
    const auto *osf_255 = buffer.data(osf + 255);
    const auto *osf_256 = buffer.data(osf + 256);
    const auto *osf_257 = buffer.data(osf + 257);
    const auto *osf_258 = buffer.data(osf + 258);
    const auto *osf_259 = buffer.data(osf + 259);
    const auto *osf_260 = buffer.data(osf + 260);
    const auto *osf_262 = buffer.data(osf + 262);
    const auto *osf_266 = buffer.data(osf + 266);
    const auto *osf_267 = buffer.data(osf + 267);
    const auto *osf_268 = buffer.data(osf + 268);
    const auto *osf_269 = buffer.data(osf + 269);
    const auto *osf_270 = buffer.data(osf + 270);
    const auto *osf_271 = buffer.data(osf + 271);
    const auto *osf_272 = buffer.data(osf + 272);
    const auto *osf_275 = buffer.data(osf + 275);
    const auto *osf_276 = buffer.data(osf + 276);
    const auto *osf_277 = buffer.data(osf + 277);
    const auto *osf_278 = buffer.data(osf + 278);
    const auto *osf_279 = buffer.data(osf + 279);
    const auto *osf_280 = buffer.data(osf + 280);
    const auto *osf_281 = buffer.data(osf + 281);
    const auto *osf_282 = buffer.data(osf + 282);
    const auto *osf_283 = buffer.data(osf + 283);
    const auto *osf_286 = buffer.data(osf + 286);
    const auto *osf_287 = buffer.data(osf + 287);
    const auto *osf_288 = buffer.data(osf + 288);
    const auto *osf_289 = buffer.data(osf + 289);
    const auto *osf_290 = buffer.data(osf + 290);
    const auto *osf_292 = buffer.data(osf + 292);
    const auto *osf_295 = buffer.data(osf + 295);
    const auto *osf_296 = buffer.data(osf + 296);
    const auto *osf_297 = buffer.data(osf + 297);
    const auto *osf_298 = buffer.data(osf + 298);
    const auto *osf_299 = buffer.data(osf + 299);
    const auto *osf_300 = buffer.data(osf + 300);
    const auto *osf_302 = buffer.data(osf + 302);
    const auto *osf_303 = buffer.data(osf + 303);
    const auto *osf_305 = buffer.data(osf + 305);
    const auto *osf_306 = buffer.data(osf + 306);
    const auto *osf_307 = buffer.data(osf + 307);
    const auto *osf_308 = buffer.data(osf + 308);
    const auto *osf_309 = buffer.data(osf + 309);
    const auto *osf_310 = buffer.data(osf + 310);
    const auto *osf_312 = buffer.data(osf + 312);
    const auto *osf_313 = buffer.data(osf + 313);
    const auto *osf_315 = buffer.data(osf + 315);
    const auto *osf_316 = buffer.data(osf + 316);
    const auto *osf_317 = buffer.data(osf + 317);
    const auto *osf_318 = buffer.data(osf + 318);
    const auto *osf_319 = buffer.data(osf + 319);
    const auto *osf_320 = buffer.data(osf + 320);
    const auto *osf_322 = buffer.data(osf + 322);
    const auto *osf_323 = buffer.data(osf + 323);

#pragma omp simd aligned(t_370, t_371, t_372, pc_y, pc_z, nsf_176, nsf_186, nsf_188, osd0_147, \
                         osd0_149, osd1_147, osd1_149, osf_246, \
                         osf_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = f_14 * nsf_186[k]
                   + f_1 * osd0_147[k]
                   - f_2 * osd1_147[k]
                   + f_3 * pc_y[k] * osf_246[k];

        t_371[k] = f_14 * nsf_176[k]
                   + f_3 * pc_z[k] * osf_246[k];

        t_372[k] = f_14 * nsf_188[k]
                   + f_4 * osd0_149[k]
                   - f_5 * osd1_149[k]
                   + f_3 * pc_y[k] * osf_248[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, pc_x, pc_y, pc_z, nsf_179, nsf_189, nsf_250, \
                         osd0_149, osd0_150, osd1_149, osd1_150, osf_249, \
                         osf_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_14 * nsf_189[k]
                   + f_3 * pc_y[k] * osf_249[k];

        t_374[k] = f_14 * nsf_179[k]
                   + f_1 * osd0_149[k]
                   - f_2 * osd1_149[k]
                   + f_3 * pc_z[k] * osf_249[k];

        t_375[k] = f_18 * nsf_250[k]
                   + f_1 * osd0_150[k]
                   - f_2 * osd1_150[k]
                   + f_3 * pc_x[k] * osf_250[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, t_379, pc_x, pc_y, pc_z, nsf_180, nsf_190, \
                         nsf_192, nsf_253, osd0_153, osd1_153, osf_250, osf_252, \
                         osf_253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_8 * nsf_190[k]
                   + f_3 * pc_y[k] * osf_250[k];

        t_377[k] = f_16 * nsf_180[k]
                   + f_3 * pc_z[k] * osf_250[k];

        t_378[k] = f_18 * nsf_253[k]
                   + f_4 * osd0_153[k]
                   - f_5 * osd1_153[k]
                   + f_3 * pc_x[k] * osf_253[k];

        t_379[k] = f_8 * nsf_192[k]
                   + f_3 * pc_y[k] * osf_252[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, pc_x, nsf_255, nsf_256, nsf_257, nsf_258, \
                         osd0_155, osd1_155, osf_255, osf_256, osf_257, \
                         osf_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = f_18 * nsf_255[k]
                   + f_4 * osd0_155[k]
                   - f_5 * osd1_155[k]
                   + f_3 * pc_x[k] * osf_255[k];

        t_381[k] = f_18 * nsf_256[k]
                   + f_3 * pc_x[k] * osf_256[k];

        t_382[k] = f_18 * nsf_257[k]
                   + f_3 * pc_x[k] * osf_257[k];

        t_383[k] = f_18 * nsf_258[k]
                   + f_3 * pc_x[k] * osf_258[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, pc_x, pc_y, pc_z, nsf_186, nsf_196, nsf_259, \
                         osd0_153, osd1_153, osf_256, osf_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = f_18 * nsf_259[k]
                   + f_3 * pc_x[k] * osf_259[k];

        t_385[k] = f_8 * nsf_196[k]
                   + f_1 * osd0_153[k]
                   - f_2 * osd1_153[k]
                   + f_3 * pc_y[k] * osf_256[k];

        t_386[k] = f_16 * nsf_186[k]
                   + f_3 * pc_z[k] * osf_256[k];
    }

#pragma omp simd aligned(t_387, t_388, t_389, t_390, pa_y, pc_y, pc_z, nsg0_300, nsf_189, \
                         nsf_198, nsf_199, nsg1_300, osd0_155, osd1_155, osf_258, \
                         osf_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_387[k] = f_8 * nsf_198[k]
                   + f_4 * osd0_155[k]
                   - f_5 * osd1_155[k]
                   + f_3 * pc_y[k] * osf_258[k];

        t_388[k] = f_8 * nsf_199[k]
                   + f_3 * pc_y[k] * osf_259[k];

        t_389[k] = f_16 * nsf_189[k]
                   + f_1 * osd0_155[k]
                   - f_2 * osd1_155[k]
                   + f_3 * pc_z[k] * osf_259[k];

        t_390[k] = pa_y[k] * nsg0_300[k]
                   - f_6 * pc_y[k] * nsg1_300[k];
    }

#pragma omp simd aligned(t_391, t_392, t_393, t_394, pa_y, pc_y, pc_z, nsg0_303, nsf_190, \
                         nsf_200, nsf_201, nsf_202, nsg1_303, osf_260, \
                         osf_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_391[k] = f_7 * nsf_200[k]
                   + f_3 * pc_y[k] * osf_260[k];

        t_392[k] = f_18 * nsf_190[k]
                   + f_3 * pc_z[k] * osf_260[k];

        t_393[k] = pa_y[k] * nsg0_303[k]
                   + f_8 * nsf_201[k]
                   - f_6 * pc_y[k] * nsg1_303[k];

        t_394[k] = f_7 * nsf_202[k]
                   + f_3 * pc_y[k] * osf_262[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, pa_y, pc_x, pc_y, nsg0_305, nsf_266, \
                         nsf_267, nsf_268, nsg1_305, osf_266, osf_267, \
                         osf_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = pa_y[k] * nsg0_305[k]
                   - f_6 * pc_y[k] * nsg1_305[k];

        t_396[k] = f_18 * nsf_266[k]
                   + f_3 * pc_x[k] * osf_266[k];

        t_397[k] = f_18 * nsf_267[k]
                   + f_3 * pc_x[k] * osf_267[k];

        t_398[k] = f_18 * nsf_268[k]
                   + f_3 * pc_x[k] * osf_268[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, pc_x, pc_y, pc_z, nsf_196, nsf_206, nsf_269, \
                         osd0_159, osd1_159, osf_266, osf_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_399[k] = f_18 * nsf_269[k]
                   + f_3 * pc_x[k] * osf_269[k];

        t_400[k] = f_7 * nsf_206[k]
                   + f_1 * osd0_159[k]
                   - f_2 * osd1_159[k]
                   + f_3 * pc_y[k] * osf_266[k];

        t_401[k] = f_18 * nsf_196[k]
                   + f_3 * pc_z[k] * osf_266[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, pa_y, pc_y, nsg0_314, nsf_208, nsf_209, \
                         nsg1_314, osd0_161, osd1_161, osf_268, \
                         osf_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = f_7 * nsf_208[k]
                   + f_4 * osd0_161[k]
                   - f_5 * osd1_161[k]
                   + f_3 * pc_y[k] * osf_268[k];

        t_403[k] = f_7 * nsf_209[k]
                   + f_3 * pc_y[k] * osf_269[k];

        t_404[k] = pa_y[k] * nsg0_314[k]
                   - f_6 * pc_y[k] * nsg1_314[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, pc_x, pc_y, pc_z, nsf_200, \
                         nsf_270, osd0_162, osd1_162, osf_270, osf_271, \
                         osf_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = f_18 * nsf_270[k]
                   + f_1 * osd0_162[k]
                   - f_2 * osd1_162[k]
                   + f_3 * pc_x[k] * osf_270[k];

        t_406[k] = f_3 * pc_y[k] * osf_270[k];

        t_407[k] = f_17 * nsf_200[k]
                   + f_3 * pc_z[k] * osf_270[k];

        t_408[k] = f_4 * osd0_162[k]
                   - f_5 * osd1_162[k]
                   + f_3 * pc_y[k] * osf_271[k];

        t_409[k] = f_3 * pc_y[k] * osf_272[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, pc_x, pc_y, nsf_275, nsf_276, nsf_277, \
                         osd0_167, osd1_167, osf_275, osf_276, \
                         osf_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_18 * nsf_275[k]
                   + f_4 * osd0_167[k]
                   - f_5 * osd1_167[k]
                   + f_3 * pc_x[k] * osf_275[k];

        t_411[k] = f_18 * nsf_276[k]
                   + f_3 * pc_x[k] * osf_276[k];

        t_412[k] = f_18 * nsf_277[k]
                   + f_3 * pc_x[k] * osf_277[k];

        t_413[k] = f_3 * pc_y[k] * osf_275[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pc_x, pc_y, nsf_279, osd0_165, osd0_166, \
                         osd1_165, osd1_166, osf_276, osf_277, \
                         osf_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_18 * nsf_279[k]
                   + f_3 * pc_x[k] * osf_279[k];

        t_415[k] = f_1 * osd0_165[k]
                   - f_2 * osd1_165[k]
                   + f_3 * pc_y[k] * osf_276[k];

        t_416[k] = f_10 * osd0_166[k]
                   - f_11 * osd1_166[k]
                   + f_3 * pc_y[k] * osf_277[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, t_420, pc_x, pc_y, pc_z, nsf_209, nsf_280, \
                         osd0_167, osd0_168, osd1_167, osd1_168, osf_278, osf_279, \
                         osf_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_4 * osd0_167[k]
                   - f_5 * osd1_167[k]
                   + f_3 * pc_y[k] * osf_278[k];

        t_418[k] = f_3 * pc_y[k] * osf_279[k];

        t_419[k] = f_17 * nsf_209[k]
                   + f_1 * osd0_167[k]
                   - f_2 * osd1_167[k]
                   + f_3 * pc_z[k] * osf_279[k];

        t_420[k] = f_16 * nsf_280[k]
                   + f_1 * osd0_168[k]
                   - f_2 * osd1_168[k]
                   + f_3 * pc_x[k] * osf_280[k];
    }

#pragma omp simd aligned(t_421, t_422, t_423, t_424, pc_x, pc_y, pc_z, nsf_210, nsf_283, \
                         osd0_171, osd1_171, osf_280, osf_281, \
                         osf_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_421[k] = f_15 * nsf_210[k]
                   + f_3 * pc_y[k] * osf_280[k];

        t_422[k] = f_3 * pc_z[k] * osf_280[k];

        t_423[k] = f_16 * nsf_283[k]
                   + f_4 * osd0_171[k]
                   - f_5 * osd1_171[k]
                   + f_3 * pc_x[k] * osf_283[k];

        t_424[k] = f_3 * pc_z[k] * osf_281[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, t_428, pc_x, pc_z, nsf_286, nsf_288, osd0_168, \
                         osd1_168, osf_282, osf_283, osf_286, osf_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = f_4 * osd0_168[k]
                   - f_5 * osd1_168[k]
                   + f_3 * pc_z[k] * osf_282[k];

        t_426[k] = f_16 * nsf_286[k]
                   + f_3 * pc_x[k] * osf_286[k];

        t_427[k] = f_3 * pc_z[k] * osf_283[k];

        t_428[k] = f_16 * nsf_288[k]
                   + f_3 * pc_x[k] * osf_288[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, t_432, t_433, pc_x, pc_y, pc_z, nsf_216, \
                         nsf_219, nsf_289, osd0_171, osd1_171, osf_286, osf_287, \
                         osf_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_16 * nsf_289[k]
                   + f_3 * pc_x[k] * osf_289[k];

        t_430[k] = f_15 * nsf_216[k]
                   + f_1 * osd0_171[k]
                   - f_2 * osd1_171[k]
                   + f_3 * pc_y[k] * osf_286[k];

        t_431[k] = f_3 * pc_z[k] * osf_286[k];

        t_432[k] = f_4 * osd0_171[k]
                   - f_5 * osd1_171[k]
                   + f_3 * pc_z[k] * osf_287[k];

        t_433[k] = f_15 * nsf_219[k]
                   + f_3 * pc_y[k] * osf_289[k];
    }

#pragma omp simd aligned(t_434, t_435, t_436, t_437, pa_z, pc_y, pc_z, nsg0_315, nsf_210, \
                         nsf_220, nsg1_315, osd0_173, osd1_173, osf_289, \
                         osf_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_434[k] = f_1 * osd0_173[k]
                   - f_2 * osd1_173[k]
                   + f_3 * pc_z[k] * osf_289[k];

        t_435[k] = pa_z[k] * nsg0_315[k]
                   - f_6 * pc_z[k] * nsg1_315[k];

        t_436[k] = f_17 * nsf_220[k]
                   + f_3 * pc_y[k] * osf_290[k];

        t_437[k] = f_7 * nsf_210[k]
                   + f_3 * pc_z[k] * osf_290[k];
    }

#pragma omp simd aligned(t_438, t_439, t_440, pa_z, pc_x, pc_y, pc_z, nsg0_318, nsf_222, \
                         nsf_295, nsg1_318, osd0_179, osd1_179, osf_292, \
                         osf_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_438[k] = pa_z[k] * nsg0_318[k]
                   - f_6 * pc_z[k] * nsg1_318[k];

        t_439[k] = f_17 * nsf_222[k]
                   + f_3 * pc_y[k] * osf_292[k];

        t_440[k] = f_16 * nsf_295[k]
                   + f_4 * osd0_179[k]
                   - f_5 * osd1_179[k]
                   + f_3 * pc_x[k] * osf_295[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, pc_x, nsf_296, nsf_297, nsf_298, nsf_299, \
                         osf_296, osf_297, osf_298, osf_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_16 * nsf_296[k]
                   + f_3 * pc_x[k] * osf_296[k];

        t_442[k] = f_16 * nsf_297[k]
                   + f_3 * pc_x[k] * osf_297[k];

        t_443[k] = f_16 * nsf_298[k]
                   + f_3 * pc_x[k] * osf_298[k];

        t_444[k] = f_16 * nsf_299[k]
                   + f_3 * pc_x[k] * osf_299[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, pa_z, pc_y, pc_z, nsg0_325, nsf_216, nsf_228, \
                         nsg1_325, osd0_179, osd1_179, osf_296, \
                         osf_298 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = pa_z[k] * nsg0_325[k]
                   - f_6 * pc_z[k] * nsg1_325[k];

        t_446[k] = f_7 * nsf_216[k]
                   + f_3 * pc_z[k] * osf_296[k];

        t_447[k] = f_17 * nsf_228[k]
                   + f_4 * osd0_179[k]
                   - f_5 * osd1_179[k]
                   + f_3 * pc_y[k] * osf_298[k];
    }

#pragma omp simd aligned(t_448, t_449, t_450, pc_x, pc_y, pc_z, nsf_219, nsf_229, nsf_300, \
                         osd0_179, osd0_180, osd1_179, osd1_180, osf_299, \
                         osf_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_448[k] = f_17 * nsf_229[k]
                   + f_3 * pc_y[k] * osf_299[k];

        t_449[k] = f_7 * nsf_219[k]
                   + f_1 * osd0_179[k]
                   - f_2 * osd1_179[k]
                   + f_3 * pc_z[k] * osf_299[k];

        t_450[k] = f_16 * nsf_300[k]
                   + f_1 * osd0_180[k]
                   - f_2 * osd1_180[k]
                   + f_3 * pc_x[k] * osf_300[k];
    }

#pragma omp simd aligned(t_451, t_452, t_453, t_454, pc_x, pc_y, pc_z, nsf_220, nsf_230, \
                         nsf_232, nsf_303, osd0_183, osd1_183, osf_300, osf_302, \
                         osf_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_451[k] = f_18 * nsf_230[k]
                   + f_3 * pc_y[k] * osf_300[k];

        t_452[k] = f_8 * nsf_220[k]
                   + f_3 * pc_z[k] * osf_300[k];

        t_453[k] = f_16 * nsf_303[k]
                   + f_4 * osd0_183[k]
                   - f_5 * osd1_183[k]
                   + f_3 * pc_x[k] * osf_303[k];

        t_454[k] = f_18 * nsf_232[k]
                   + f_3 * pc_y[k] * osf_302[k];
    }

#pragma omp simd aligned(t_455, t_456, t_457, t_458, pc_x, nsf_305, nsf_306, nsf_307, nsf_308, \
                         osd0_185, osd1_185, osf_305, osf_306, osf_307, \
                         osf_308 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_455[k] = f_16 * nsf_305[k]
                   + f_4 * osd0_185[k]
                   - f_5 * osd1_185[k]
                   + f_3 * pc_x[k] * osf_305[k];

        t_456[k] = f_16 * nsf_306[k]
                   + f_3 * pc_x[k] * osf_306[k];

        t_457[k] = f_16 * nsf_307[k]
                   + f_3 * pc_x[k] * osf_307[k];

        t_458[k] = f_16 * nsf_308[k]
                   + f_3 * pc_x[k] * osf_308[k];
    }

#pragma omp simd aligned(t_459, t_460, t_461, pc_x, pc_y, pc_z, nsf_226, nsf_236, nsf_309, \
                         osd0_183, osd1_183, osf_306, osf_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_459[k] = f_16 * nsf_309[k]
                   + f_3 * pc_x[k] * osf_309[k];

        t_460[k] = f_18 * nsf_236[k]
                   + f_1 * osd0_183[k]
                   - f_2 * osd1_183[k]
                   + f_3 * pc_y[k] * osf_306[k];

        t_461[k] = f_8 * nsf_226[k]
                   + f_3 * pc_z[k] * osf_306[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, pc_y, pc_z, nsf_229, nsf_238, nsf_239, osd0_185, \
                         osd1_185, osf_308, osf_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = f_18 * nsf_238[k]
                   + f_4 * osd0_185[k]
                   - f_5 * osd1_185[k]
                   + f_3 * pc_y[k] * osf_308[k];

        t_463[k] = f_18 * nsf_239[k]
                   + f_3 * pc_y[k] * osf_309[k];

        t_464[k] = f_8 * nsf_229[k]
                   + f_1 * osd0_185[k]
                   - f_2 * osd1_185[k]
                   + f_3 * pc_z[k] * osf_309[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, pc_x, pc_y, pc_z, nsf_230, nsf_240, nsf_310, \
                         osd0_186, osd1_186, osf_310 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = f_16 * nsf_310[k]
                   + f_1 * osd0_186[k]
                   - f_2 * osd1_186[k]
                   + f_3 * pc_x[k] * osf_310[k];

        t_466[k] = f_16 * nsf_240[k]
                   + f_3 * pc_y[k] * osf_310[k];

        t_467[k] = f_14 * nsf_230[k]
                   + f_3 * pc_z[k] * osf_310[k];
    }

#pragma omp simd aligned(t_468, t_469, t_470, pc_x, pc_y, nsf_242, nsf_313, nsf_315, osd0_189, \
                         osd0_191, osd1_189, osd1_191, osf_312, osf_313, \
                         osf_315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_468[k] = f_16 * nsf_313[k]
                   + f_4 * osd0_189[k]
                   - f_5 * osd1_189[k]
                   + f_3 * pc_x[k] * osf_313[k];

        t_469[k] = f_16 * nsf_242[k]
                   + f_3 * pc_y[k] * osf_312[k];

        t_470[k] = f_16 * nsf_315[k]
                   + f_4 * osd0_191[k]
                   - f_5 * osd1_191[k]
                   + f_3 * pc_x[k] * osf_315[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, t_474, pc_x, nsf_316, nsf_317, nsf_318, nsf_319, \
                         osf_316, osf_317, osf_318, osf_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = f_16 * nsf_316[k]
                   + f_3 * pc_x[k] * osf_316[k];

        t_472[k] = f_16 * nsf_317[k]
                   + f_3 * pc_x[k] * osf_317[k];

        t_473[k] = f_16 * nsf_318[k]
                   + f_3 * pc_x[k] * osf_318[k];

        t_474[k] = f_16 * nsf_319[k]
                   + f_3 * pc_x[k] * osf_319[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, pc_y, pc_z, nsf_236, nsf_246, nsf_248, osd0_189, \
                         osd0_191, osd1_189, osd1_191, osf_316, \
                         osf_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = f_16 * nsf_246[k]
                   + f_1 * osd0_189[k]
                   - f_2 * osd1_189[k]
                   + f_3 * pc_y[k] * osf_316[k];

        t_476[k] = f_14 * nsf_236[k]
                   + f_3 * pc_z[k] * osf_316[k];

        t_477[k] = f_16 * nsf_248[k]
                   + f_4 * osd0_191[k]
                   - f_5 * osd1_191[k]
                   + f_3 * pc_y[k] * osf_318[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, pc_x, pc_y, pc_z, nsf_239, nsf_249, nsf_320, \
                         osd0_191, osd0_192, osd1_191, osd1_192, osf_319, \
                         osf_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = f_16 * nsf_249[k]
                   + f_3 * pc_y[k] * osf_319[k];

        t_479[k] = f_14 * nsf_239[k]
                   + f_1 * osd0_191[k]
                   - f_2 * osd1_191[k]
                   + f_3 * pc_z[k] * osf_319[k];

        t_480[k] = f_16 * nsf_320[k]
                   + f_1 * osd0_192[k]
                   - f_2 * osd1_192[k]
                   + f_3 * pc_x[k] * osf_320[k];
    }

#pragma omp simd aligned(t_481, t_482, t_483, t_484, pc_x, pc_y, pc_z, nsf_240, nsf_250, \
                         nsf_252, nsf_323, osd0_195, osd1_195, osf_320, osf_322, \
                         osf_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_481[k] = f_14 * nsf_250[k]
                   + f_3 * pc_y[k] * osf_320[k];

        t_482[k] = f_16 * nsf_240[k]
                   + f_3 * pc_z[k] * osf_320[k];

        t_483[k] = f_16 * nsf_323[k]
                   + f_4 * osd0_195[k]
                   - f_5 * osd1_195[k]
                   + f_3 * pc_x[k] * osf_323[k];

        t_484[k] = f_14 * nsf_252[k]
                   + f_3 * pc_y[k] * osf_322[k];
    }
}

static auto
compute_prim_osg_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsg0,
                                                          const size_t nsf, const size_t nsg1,
                                                          const size_t osd0, const size_t osd1,
                                                          const size_t osf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_13 = 4.0 / q;
    const auto f_14 = 1.5 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.5 / q;

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
    auto *t_593 = buffer.data(target + 593);
    auto *t_594 = buffer.data(target + 594);
    auto *t_595 = buffer.data(target + 595);
    auto *t_596 = buffer.data(target + 596);
    auto *t_597 = buffer.data(target + 597);
    auto *t_598 = buffer.data(target + 598);
    auto *t_599 = buffer.data(target + 599);
    auto *t_600 = buffer.data(target + 600);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsg0_405 = buffer.data(nsg0 + 405);
    const auto *nsg0_408 = buffer.data(nsg0 + 408);
    const auto *nsg0_410 = buffer.data(nsg0 + 410);
    const auto *nsg0_419 = buffer.data(nsg0 + 419);
    const auto *nsg0_420 = buffer.data(nsg0 + 420);
    const auto *nsg0_423 = buffer.data(nsg0 + 423);
    const auto *nsg0_430 = buffer.data(nsg0 + 430);

    const auto *nsf_246 = buffer.data(nsf + 246);
    const auto *nsf_249 = buffer.data(nsf + 249);
    const auto *nsf_250 = buffer.data(nsf + 250);
    const auto *nsf_256 = buffer.data(nsf + 256);
    const auto *nsf_258 = buffer.data(nsf + 258);
    const auto *nsf_259 = buffer.data(nsf + 259);
    const auto *nsf_260 = buffer.data(nsf + 260);
    const auto *nsf_262 = buffer.data(nsf + 262);
    const auto *nsf_266 = buffer.data(nsf + 266);
    const auto *nsf_268 = buffer.data(nsf + 268);
    const auto *nsf_269 = buffer.data(nsf + 269);
    const auto *nsf_270 = buffer.data(nsf + 270);
    const auto *nsf_271 = buffer.data(nsf + 271);
    const auto *nsf_272 = buffer.data(nsf + 272);
    const auto *nsf_276 = buffer.data(nsf + 276);
    const auto *nsf_278 = buffer.data(nsf + 278);
    const auto *nsf_279 = buffer.data(nsf + 279);
    const auto *nsf_280 = buffer.data(nsf + 280);
    const auto *nsf_286 = buffer.data(nsf + 286);
    const auto *nsf_289 = buffer.data(nsf + 289);
    const auto *nsf_290 = buffer.data(nsf + 290);
    const auto *nsf_292 = buffer.data(nsf + 292);
    const auto *nsf_296 = buffer.data(nsf + 296);
    const auto *nsf_298 = buffer.data(nsf + 298);
    const auto *nsf_299 = buffer.data(nsf + 299);
    const auto *nsf_300 = buffer.data(nsf + 300);
    const auto *nsf_302 = buffer.data(nsf + 302);
    const auto *nsf_306 = buffer.data(nsf + 306);
    const auto *nsf_308 = buffer.data(nsf + 308);
    const auto *nsf_309 = buffer.data(nsf + 309);
    const auto *nsf_310 = buffer.data(nsf + 310);
    const auto *nsf_312 = buffer.data(nsf + 312);
    const auto *nsf_316 = buffer.data(nsf + 316);
    const auto *nsf_318 = buffer.data(nsf + 318);
    const auto *nsf_319 = buffer.data(nsf + 319);
    const auto *nsf_325 = buffer.data(nsf + 325);
    const auto *nsf_326 = buffer.data(nsf + 326);
    const auto *nsf_327 = buffer.data(nsf + 327);
    const auto *nsf_328 = buffer.data(nsf + 328);
    const auto *nsf_329 = buffer.data(nsf + 329);
    const auto *nsf_330 = buffer.data(nsf + 330);
    const auto *nsf_333 = buffer.data(nsf + 333);
    const auto *nsf_335 = buffer.data(nsf + 335);
    const auto *nsf_336 = buffer.data(nsf + 336);
    const auto *nsf_337 = buffer.data(nsf + 337);
    const auto *nsf_338 = buffer.data(nsf + 338);
    const auto *nsf_339 = buffer.data(nsf + 339);
    const auto *nsf_346 = buffer.data(nsf + 346);
    const auto *nsf_347 = buffer.data(nsf + 347);
    const auto *nsf_348 = buffer.data(nsf + 348);
    const auto *nsf_349 = buffer.data(nsf + 349);
    const auto *nsf_350 = buffer.data(nsf + 350);
    const auto *nsf_355 = buffer.data(nsf + 355);
    const auto *nsf_356 = buffer.data(nsf + 356);
    const auto *nsf_357 = buffer.data(nsf + 357);
    const auto *nsf_359 = buffer.data(nsf + 359);
    const auto *nsf_360 = buffer.data(nsf + 360);
    const auto *nsf_363 = buffer.data(nsf + 363);
    const auto *nsf_366 = buffer.data(nsf + 366);
    const auto *nsf_368 = buffer.data(nsf + 368);
    const auto *nsf_369 = buffer.data(nsf + 369);
    const auto *nsf_375 = buffer.data(nsf + 375);
    const auto *nsf_376 = buffer.data(nsf + 376);
    const auto *nsf_377 = buffer.data(nsf + 377);
    const auto *nsf_378 = buffer.data(nsf + 378);
    const auto *nsf_379 = buffer.data(nsf + 379);
    const auto *nsf_380 = buffer.data(nsf + 380);
    const auto *nsf_383 = buffer.data(nsf + 383);
    const auto *nsf_385 = buffer.data(nsf + 385);
    const auto *nsf_386 = buffer.data(nsf + 386);
    const auto *nsf_387 = buffer.data(nsf + 387);
    const auto *nsf_388 = buffer.data(nsf + 388);
    const auto *nsf_389 = buffer.data(nsf + 389);
    const auto *nsf_390 = buffer.data(nsf + 390);
    const auto *nsf_393 = buffer.data(nsf + 393);
    const auto *nsf_395 = buffer.data(nsf + 395);
    const auto *nsf_396 = buffer.data(nsf + 396);
    const auto *nsf_397 = buffer.data(nsf + 397);
    const auto *nsf_398 = buffer.data(nsf + 398);
    const auto *nsf_399 = buffer.data(nsf + 399);
    const auto *nsf_400 = buffer.data(nsf + 400);

    const auto *nsg1_405 = buffer.data(nsg1 + 405);
    const auto *nsg1_408 = buffer.data(nsg1 + 408);
    const auto *nsg1_410 = buffer.data(nsg1 + 410);
    const auto *nsg1_419 = buffer.data(nsg1 + 419);
    const auto *nsg1_420 = buffer.data(nsg1 + 420);
    const auto *nsg1_423 = buffer.data(nsg1 + 423);
    const auto *nsg1_430 = buffer.data(nsg1 + 430);

    const auto *osd0_195 = buffer.data(osd0 + 195);
    const auto *osd0_197 = buffer.data(osd0 + 197);
    const auto *osd0_198 = buffer.data(osd0 + 198);
    const auto *osd0_201 = buffer.data(osd0 + 201);
    const auto *osd0_203 = buffer.data(osd0 + 203);
    const auto *osd0_207 = buffer.data(osd0 + 207);
    const auto *osd0_209 = buffer.data(osd0 + 209);
    const auto *osd0_210 = buffer.data(osd0 + 210);
    const auto *osd0_213 = buffer.data(osd0 + 213);
    const auto *osd0_214 = buffer.data(osd0 + 214);
    const auto *osd0_215 = buffer.data(osd0 + 215);
    const auto *osd0_216 = buffer.data(osd0 + 216);
    const auto *osd0_219 = buffer.data(osd0 + 219);
    const auto *osd0_221 = buffer.data(osd0 + 221);
    const auto *osd0_227 = buffer.data(osd0 + 227);
    const auto *osd0_228 = buffer.data(osd0 + 228);
    const auto *osd0_231 = buffer.data(osd0 + 231);
    const auto *osd0_233 = buffer.data(osd0 + 233);
    const auto *osd0_234 = buffer.data(osd0 + 234);
    const auto *osd0_237 = buffer.data(osd0 + 237);
    const auto *osd0_239 = buffer.data(osd0 + 239);
    const auto *osd0_240 = buffer.data(osd0 + 240);

    const auto *osd1_195 = buffer.data(osd1 + 195);
    const auto *osd1_197 = buffer.data(osd1 + 197);
    const auto *osd1_198 = buffer.data(osd1 + 198);
    const auto *osd1_201 = buffer.data(osd1 + 201);
    const auto *osd1_203 = buffer.data(osd1 + 203);
    const auto *osd1_207 = buffer.data(osd1 + 207);
    const auto *osd1_209 = buffer.data(osd1 + 209);
    const auto *osd1_210 = buffer.data(osd1 + 210);
    const auto *osd1_213 = buffer.data(osd1 + 213);
    const auto *osd1_214 = buffer.data(osd1 + 214);
    const auto *osd1_215 = buffer.data(osd1 + 215);
    const auto *osd1_216 = buffer.data(osd1 + 216);
    const auto *osd1_219 = buffer.data(osd1 + 219);
    const auto *osd1_221 = buffer.data(osd1 + 221);
    const auto *osd1_227 = buffer.data(osd1 + 227);
    const auto *osd1_228 = buffer.data(osd1 + 228);
    const auto *osd1_231 = buffer.data(osd1 + 231);
    const auto *osd1_233 = buffer.data(osd1 + 233);
    const auto *osd1_234 = buffer.data(osd1 + 234);
    const auto *osd1_237 = buffer.data(osd1 + 237);
    const auto *osd1_239 = buffer.data(osd1 + 239);
    const auto *osd1_240 = buffer.data(osd1 + 240);

    const auto *osf_325 = buffer.data(osf + 325);
    const auto *osf_326 = buffer.data(osf + 326);
    const auto *osf_327 = buffer.data(osf + 327);
    const auto *osf_328 = buffer.data(osf + 328);
    const auto *osf_329 = buffer.data(osf + 329);
    const auto *osf_330 = buffer.data(osf + 330);
    const auto *osf_332 = buffer.data(osf + 332);
    const auto *osf_333 = buffer.data(osf + 333);
    const auto *osf_335 = buffer.data(osf + 335);
    const auto *osf_336 = buffer.data(osf + 336);
    const auto *osf_337 = buffer.data(osf + 337);
    const auto *osf_338 = buffer.data(osf + 338);
    const auto *osf_339 = buffer.data(osf + 339);
    const auto *osf_340 = buffer.data(osf + 340);
    const auto *osf_342 = buffer.data(osf + 342);
    const auto *osf_346 = buffer.data(osf + 346);
    const auto *osf_347 = buffer.data(osf + 347);
    const auto *osf_348 = buffer.data(osf + 348);
    const auto *osf_349 = buffer.data(osf + 349);
    const auto *osf_350 = buffer.data(osf + 350);
    const auto *osf_351 = buffer.data(osf + 351);
    const auto *osf_352 = buffer.data(osf + 352);
    const auto *osf_355 = buffer.data(osf + 355);
    const auto *osf_356 = buffer.data(osf + 356);
    const auto *osf_357 = buffer.data(osf + 357);
    const auto *osf_358 = buffer.data(osf + 358);
    const auto *osf_359 = buffer.data(osf + 359);
    const auto *osf_360 = buffer.data(osf + 360);
    const auto *osf_361 = buffer.data(osf + 361);
    const auto *osf_362 = buffer.data(osf + 362);
    const auto *osf_363 = buffer.data(osf + 363);
    const auto *osf_366 = buffer.data(osf + 366);
    const auto *osf_367 = buffer.data(osf + 367);
    const auto *osf_368 = buffer.data(osf + 368);
    const auto *osf_369 = buffer.data(osf + 369);
    const auto *osf_370 = buffer.data(osf + 370);
    const auto *osf_372 = buffer.data(osf + 372);
    const auto *osf_375 = buffer.data(osf + 375);
    const auto *osf_376 = buffer.data(osf + 376);
    const auto *osf_377 = buffer.data(osf + 377);
    const auto *osf_378 = buffer.data(osf + 378);
    const auto *osf_379 = buffer.data(osf + 379);
    const auto *osf_380 = buffer.data(osf + 380);
    const auto *osf_382 = buffer.data(osf + 382);
    const auto *osf_383 = buffer.data(osf + 383);
    const auto *osf_385 = buffer.data(osf + 385);
    const auto *osf_386 = buffer.data(osf + 386);
    const auto *osf_387 = buffer.data(osf + 387);
    const auto *osf_388 = buffer.data(osf + 388);
    const auto *osf_389 = buffer.data(osf + 389);
    const auto *osf_390 = buffer.data(osf + 390);
    const auto *osf_392 = buffer.data(osf + 392);
    const auto *osf_393 = buffer.data(osf + 393);
    const auto *osf_395 = buffer.data(osf + 395);
    const auto *osf_396 = buffer.data(osf + 396);
    const auto *osf_397 = buffer.data(osf + 397);
    const auto *osf_398 = buffer.data(osf + 398);
    const auto *osf_399 = buffer.data(osf + 399);
    const auto *osf_400 = buffer.data(osf + 400);

#pragma omp simd aligned(t_485, t_486, t_487, t_488, pc_x, nsf_325, nsf_326, nsf_327, nsf_328, \
                         osd0_197, osd1_197, osf_325, osf_326, osf_327, \
                         osf_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = f_16 * nsf_325[k]
                   + f_4 * osd0_197[k]
                   - f_5 * osd1_197[k]
                   + f_3 * pc_x[k] * osf_325[k];

        t_486[k] = f_16 * nsf_326[k]
                   + f_3 * pc_x[k] * osf_326[k];

        t_487[k] = f_16 * nsf_327[k]
                   + f_3 * pc_x[k] * osf_327[k];

        t_488[k] = f_16 * nsf_328[k]
                   + f_3 * pc_x[k] * osf_328[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, pc_x, pc_y, pc_z, nsf_246, nsf_256, nsf_329, \
                         osd0_195, osd1_195, osf_326, osf_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = f_16 * nsf_329[k]
                   + f_3 * pc_x[k] * osf_329[k];

        t_490[k] = f_14 * nsf_256[k]
                   + f_1 * osd0_195[k]
                   - f_2 * osd1_195[k]
                   + f_3 * pc_y[k] * osf_326[k];

        t_491[k] = f_16 * nsf_246[k]
                   + f_3 * pc_z[k] * osf_326[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, pc_y, pc_z, nsf_249, nsf_258, nsf_259, osd0_197, \
                         osd1_197, osf_328, osf_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = f_14 * nsf_258[k]
                   + f_4 * osd0_197[k]
                   - f_5 * osd1_197[k]
                   + f_3 * pc_y[k] * osf_328[k];

        t_493[k] = f_14 * nsf_259[k]
                   + f_3 * pc_y[k] * osf_329[k];

        t_494[k] = f_16 * nsf_249[k]
                   + f_1 * osd0_197[k]
                   - f_2 * osd1_197[k]
                   + f_3 * pc_z[k] * osf_329[k];
    }

#pragma omp simd aligned(t_495, t_496, t_497, pc_x, pc_y, pc_z, nsf_250, nsf_260, nsf_330, \
                         osd0_198, osd1_198, osf_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = f_16 * nsf_330[k]
                   + f_1 * osd0_198[k]
                   - f_2 * osd1_198[k]
                   + f_3 * pc_x[k] * osf_330[k];

        t_496[k] = f_8 * nsf_260[k]
                   + f_3 * pc_y[k] * osf_330[k];

        t_497[k] = f_18 * nsf_250[k]
                   + f_3 * pc_z[k] * osf_330[k];
    }

#pragma omp simd aligned(t_498, t_499, t_500, pc_x, pc_y, nsf_262, nsf_333, nsf_335, osd0_201, \
                         osd0_203, osd1_201, osd1_203, osf_332, osf_333, \
                         osf_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_498[k] = f_16 * nsf_333[k]
                   + f_4 * osd0_201[k]
                   - f_5 * osd1_201[k]
                   + f_3 * pc_x[k] * osf_333[k];

        t_499[k] = f_8 * nsf_262[k]
                   + f_3 * pc_y[k] * osf_332[k];

        t_500[k] = f_16 * nsf_335[k]
                   + f_4 * osd0_203[k]
                   - f_5 * osd1_203[k]
                   + f_3 * pc_x[k] * osf_335[k];
    }

#pragma omp simd aligned(t_501, t_502, t_503, t_504, pc_x, nsf_336, nsf_337, nsf_338, nsf_339, \
                         osf_336, osf_337, osf_338, osf_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_501[k] = f_16 * nsf_336[k]
                   + f_3 * pc_x[k] * osf_336[k];

        t_502[k] = f_16 * nsf_337[k]
                   + f_3 * pc_x[k] * osf_337[k];

        t_503[k] = f_16 * nsf_338[k]
                   + f_3 * pc_x[k] * osf_338[k];

        t_504[k] = f_16 * nsf_339[k]
                   + f_3 * pc_x[k] * osf_339[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, pc_y, pc_z, nsf_256, nsf_266, nsf_268, osd0_201, \
                         osd0_203, osd1_201, osd1_203, osf_336, \
                         osf_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = f_8 * nsf_266[k]
                   + f_1 * osd0_201[k]
                   - f_2 * osd1_201[k]
                   + f_3 * pc_y[k] * osf_336[k];

        t_506[k] = f_18 * nsf_256[k]
                   + f_3 * pc_z[k] * osf_336[k];

        t_507[k] = f_8 * nsf_268[k]
                   + f_4 * osd0_203[k]
                   - f_5 * osd1_203[k]
                   + f_3 * pc_y[k] * osf_338[k];
    }

#pragma omp simd aligned(t_508, t_509, t_510, t_511, pa_y, pc_y, pc_z, nsg0_405, nsf_259, \
                         nsf_269, nsf_270, nsg1_405, osd0_203, osd1_203, osf_339, \
                         osf_340 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_508[k] = f_8 * nsf_269[k]
                   + f_3 * pc_y[k] * osf_339[k];

        t_509[k] = f_18 * nsf_259[k]
                   + f_1 * osd0_203[k]
                   - f_2 * osd1_203[k]
                   + f_3 * pc_z[k] * osf_339[k];

        t_510[k] = pa_y[k] * nsg0_405[k]
                   - f_6 * pc_y[k] * nsg1_405[k];

        t_511[k] = f_7 * nsf_270[k]
                   + f_3 * pc_y[k] * osf_340[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, t_515, pa_y, pc_y, pc_z, nsg0_408, nsg0_410, \
                         nsf_260, nsf_271, nsf_272, nsg1_408, nsg1_410, osf_340, \
                         osf_342 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = f_17 * nsf_260[k]
                   + f_3 * pc_z[k] * osf_340[k];

        t_513[k] = pa_y[k] * nsg0_408[k]
                   + f_8 * nsf_271[k]
                   - f_6 * pc_y[k] * nsg1_408[k];

        t_514[k] = f_7 * nsf_272[k]
                   + f_3 * pc_y[k] * osf_342[k];

        t_515[k] = pa_y[k] * nsg0_410[k]
                   - f_6 * pc_y[k] * nsg1_410[k];
    }

#pragma omp simd aligned(t_516, t_517, t_518, t_519, pc_x, nsf_346, nsf_347, nsf_348, nsf_349, \
                         osf_346, osf_347, osf_348, osf_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_516[k] = f_16 * nsf_346[k]
                   + f_3 * pc_x[k] * osf_346[k];

        t_517[k] = f_16 * nsf_347[k]
                   + f_3 * pc_x[k] * osf_347[k];

        t_518[k] = f_16 * nsf_348[k]
                   + f_3 * pc_x[k] * osf_348[k];

        t_519[k] = f_16 * nsf_349[k]
                   + f_3 * pc_x[k] * osf_349[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, pc_y, pc_z, nsf_266, nsf_276, nsf_278, osd0_207, \
                         osd0_209, osd1_207, osd1_209, osf_346, \
                         osf_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = f_7 * nsf_276[k]
                   + f_1 * osd0_207[k]
                   - f_2 * osd1_207[k]
                   + f_3 * pc_y[k] * osf_346[k];

        t_521[k] = f_17 * nsf_266[k]
                   + f_3 * pc_z[k] * osf_346[k];

        t_522[k] = f_7 * nsf_278[k]
                   + f_4 * osd0_209[k]
                   - f_5 * osd1_209[k]
                   + f_3 * pc_y[k] * osf_348[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, t_526, pa_y, pc_x, pc_y, nsg0_419, nsf_279, \
                         nsf_350, nsg1_419, osd0_210, osd1_210, osf_349, \
                         osf_350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = f_7 * nsf_279[k]
                   + f_3 * pc_y[k] * osf_349[k];

        t_524[k] = pa_y[k] * nsg0_419[k]
                   - f_6 * pc_y[k] * nsg1_419[k];

        t_525[k] = f_16 * nsf_350[k]
                   + f_1 * osd0_210[k]
                   - f_2 * osd1_210[k]
                   + f_3 * pc_x[k] * osf_350[k];

        t_526[k] = f_3 * pc_y[k] * osf_350[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, pc_y, pc_z, nsf_270, osd0_210, osd1_210, \
                         osf_350, osf_351, osf_352 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = f_15 * nsf_270[k]
                   + f_3 * pc_z[k] * osf_350[k];

        t_528[k] = f_4 * osd0_210[k]
                   - f_5 * osd1_210[k]
                   + f_3 * pc_y[k] * osf_351[k];

        t_529[k] = f_3 * pc_y[k] * osf_352[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, t_533, pc_x, pc_y, nsf_355, nsf_356, nsf_357, \
                         osd0_215, osd1_215, osf_355, osf_356, \
                         osf_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_16 * nsf_355[k]
                   + f_4 * osd0_215[k]
                   - f_5 * osd1_215[k]
                   + f_3 * pc_x[k] * osf_355[k];

        t_531[k] = f_16 * nsf_356[k]
                   + f_3 * pc_x[k] * osf_356[k];

        t_532[k] = f_16 * nsf_357[k]
                   + f_3 * pc_x[k] * osf_357[k];

        t_533[k] = f_3 * pc_y[k] * osf_355[k];
    }

#pragma omp simd aligned(t_534, t_535, t_536, pc_x, pc_y, nsf_359, osd0_213, osd0_214, \
                         osd1_213, osd1_214, osf_356, osf_357, \
                         osf_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_534[k] = f_16 * nsf_359[k]
                   + f_3 * pc_x[k] * osf_359[k];

        t_535[k] = f_1 * osd0_213[k]
                   - f_2 * osd1_213[k]
                   + f_3 * pc_y[k] * osf_356[k];

        t_536[k] = f_10 * osd0_214[k]
                   - f_11 * osd1_214[k]
                   + f_3 * pc_y[k] * osf_357[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pc_x, pc_y, pc_z, nsf_279, nsf_360, \
                         osd0_215, osd0_216, osd1_215, osd1_216, osf_358, osf_359, \
                         osf_360 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = f_4 * osd0_215[k]
                   - f_5 * osd1_215[k]
                   + f_3 * pc_y[k] * osf_358[k];

        t_538[k] = f_3 * pc_y[k] * osf_359[k];

        t_539[k] = f_15 * nsf_279[k]
                   + f_1 * osd0_215[k]
                   - f_2 * osd1_215[k]
                   + f_3 * pc_z[k] * osf_359[k];

        t_540[k] = f_14 * nsf_360[k]
                   + f_1 * osd0_216[k]
                   - f_2 * osd1_216[k]
                   + f_3 * pc_x[k] * osf_360[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, t_544, pc_x, pc_y, pc_z, nsf_280, nsf_363, \
                         osd0_219, osd1_219, osf_360, osf_361, \
                         osf_363 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = f_13 * nsf_280[k]
                   + f_3 * pc_y[k] * osf_360[k];

        t_542[k] = f_3 * pc_z[k] * osf_360[k];

        t_543[k] = f_14 * nsf_363[k]
                   + f_4 * osd0_219[k]
                   - f_5 * osd1_219[k]
                   + f_3 * pc_x[k] * osf_363[k];

        t_544[k] = f_3 * pc_z[k] * osf_361[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, pc_x, pc_z, nsf_366, nsf_368, osd0_216, \
                         osd1_216, osf_362, osf_363, osf_366, osf_368 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = f_4 * osd0_216[k]
                   - f_5 * osd1_216[k]
                   + f_3 * pc_z[k] * osf_362[k];

        t_546[k] = f_14 * nsf_366[k]
                   + f_3 * pc_x[k] * osf_366[k];

        t_547[k] = f_3 * pc_z[k] * osf_363[k];

        t_548[k] = f_14 * nsf_368[k]
                   + f_3 * pc_x[k] * osf_368[k];
    }

#pragma omp simd aligned(t_549, t_550, t_551, t_552, t_553, pc_x, pc_y, pc_z, nsf_286, \
                         nsf_289, nsf_369, osd0_219, osd1_219, osf_366, osf_367, \
                         osf_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_549[k] = f_14 * nsf_369[k]
                   + f_3 * pc_x[k] * osf_369[k];

        t_550[k] = f_13 * nsf_286[k]
                   + f_1 * osd0_219[k]
                   - f_2 * osd1_219[k]
                   + f_3 * pc_y[k] * osf_366[k];

        t_551[k] = f_3 * pc_z[k] * osf_366[k];

        t_552[k] = f_4 * osd0_219[k]
                   - f_5 * osd1_219[k]
                   + f_3 * pc_z[k] * osf_367[k];

        t_553[k] = f_13 * nsf_289[k]
                   + f_3 * pc_y[k] * osf_369[k];
    }

#pragma omp simd aligned(t_554, t_555, t_556, t_557, pa_z, pc_y, pc_z, nsg0_420, nsf_280, \
                         nsf_290, nsg1_420, osd0_221, osd1_221, osf_369, \
                         osf_370 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_554[k] = f_1 * osd0_221[k]
                   - f_2 * osd1_221[k]
                   + f_3 * pc_z[k] * osf_369[k];

        t_555[k] = pa_z[k] * nsg0_420[k]
                   - f_6 * pc_z[k] * nsg1_420[k];

        t_556[k] = f_15 * nsf_290[k]
                   + f_3 * pc_y[k] * osf_370[k];

        t_557[k] = f_7 * nsf_280[k]
                   + f_3 * pc_z[k] * osf_370[k];
    }

#pragma omp simd aligned(t_558, t_559, t_560, pa_z, pc_x, pc_y, pc_z, nsg0_423, nsf_292, \
                         nsf_375, nsg1_423, osd0_227, osd1_227, osf_372, \
                         osf_375 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_558[k] = pa_z[k] * nsg0_423[k]
                   - f_6 * pc_z[k] * nsg1_423[k];

        t_559[k] = f_15 * nsf_292[k]
                   + f_3 * pc_y[k] * osf_372[k];

        t_560[k] = f_14 * nsf_375[k]
                   + f_4 * osd0_227[k]
                   - f_5 * osd1_227[k]
                   + f_3 * pc_x[k] * osf_375[k];
    }

#pragma omp simd aligned(t_561, t_562, t_563, t_564, pc_x, nsf_376, nsf_377, nsf_378, nsf_379, \
                         osf_376, osf_377, osf_378, osf_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_561[k] = f_14 * nsf_376[k]
                   + f_3 * pc_x[k] * osf_376[k];

        t_562[k] = f_14 * nsf_377[k]
                   + f_3 * pc_x[k] * osf_377[k];

        t_563[k] = f_14 * nsf_378[k]
                   + f_3 * pc_x[k] * osf_378[k];

        t_564[k] = f_14 * nsf_379[k]
                   + f_3 * pc_x[k] * osf_379[k];
    }

#pragma omp simd aligned(t_565, t_566, t_567, pa_z, pc_y, pc_z, nsg0_430, nsf_286, nsf_298, \
                         nsg1_430, osd0_227, osd1_227, osf_376, \
                         osf_378 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_565[k] = pa_z[k] * nsg0_430[k]
                   - f_6 * pc_z[k] * nsg1_430[k];

        t_566[k] = f_7 * nsf_286[k]
                   + f_3 * pc_z[k] * osf_376[k];

        t_567[k] = f_15 * nsf_298[k]
                   + f_4 * osd0_227[k]
                   - f_5 * osd1_227[k]
                   + f_3 * pc_y[k] * osf_378[k];
    }

#pragma omp simd aligned(t_568, t_569, t_570, pc_x, pc_y, pc_z, nsf_289, nsf_299, nsf_380, \
                         osd0_227, osd0_228, osd1_227, osd1_228, osf_379, \
                         osf_380 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_568[k] = f_15 * nsf_299[k]
                   + f_3 * pc_y[k] * osf_379[k];

        t_569[k] = f_7 * nsf_289[k]
                   + f_1 * osd0_227[k]
                   - f_2 * osd1_227[k]
                   + f_3 * pc_z[k] * osf_379[k];

        t_570[k] = f_14 * nsf_380[k]
                   + f_1 * osd0_228[k]
                   - f_2 * osd1_228[k]
                   + f_3 * pc_x[k] * osf_380[k];
    }

#pragma omp simd aligned(t_571, t_572, t_573, t_574, pc_x, pc_y, pc_z, nsf_290, nsf_300, \
                         nsf_302, nsf_383, osd0_231, osd1_231, osf_380, osf_382, \
                         osf_383 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_571[k] = f_17 * nsf_300[k]
                   + f_3 * pc_y[k] * osf_380[k];

        t_572[k] = f_8 * nsf_290[k]
                   + f_3 * pc_z[k] * osf_380[k];

        t_573[k] = f_14 * nsf_383[k]
                   + f_4 * osd0_231[k]
                   - f_5 * osd1_231[k]
                   + f_3 * pc_x[k] * osf_383[k];

        t_574[k] = f_17 * nsf_302[k]
                   + f_3 * pc_y[k] * osf_382[k];
    }

#pragma omp simd aligned(t_575, t_576, t_577, t_578, pc_x, nsf_385, nsf_386, nsf_387, nsf_388, \
                         osd0_233, osd1_233, osf_385, osf_386, osf_387, \
                         osf_388 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_575[k] = f_14 * nsf_385[k]
                   + f_4 * osd0_233[k]
                   - f_5 * osd1_233[k]
                   + f_3 * pc_x[k] * osf_385[k];

        t_576[k] = f_14 * nsf_386[k]
                   + f_3 * pc_x[k] * osf_386[k];

        t_577[k] = f_14 * nsf_387[k]
                   + f_3 * pc_x[k] * osf_387[k];

        t_578[k] = f_14 * nsf_388[k]
                   + f_3 * pc_x[k] * osf_388[k];
    }

#pragma omp simd aligned(t_579, t_580, t_581, pc_x, pc_y, pc_z, nsf_296, nsf_306, nsf_389, \
                         osd0_231, osd1_231, osf_386, osf_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_579[k] = f_14 * nsf_389[k]
                   + f_3 * pc_x[k] * osf_389[k];

        t_580[k] = f_17 * nsf_306[k]
                   + f_1 * osd0_231[k]
                   - f_2 * osd1_231[k]
                   + f_3 * pc_y[k] * osf_386[k];

        t_581[k] = f_8 * nsf_296[k]
                   + f_3 * pc_z[k] * osf_386[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, pc_y, pc_z, nsf_299, nsf_308, nsf_309, osd0_233, \
                         osd1_233, osf_388, osf_389 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = f_17 * nsf_308[k]
                   + f_4 * osd0_233[k]
                   - f_5 * osd1_233[k]
                   + f_3 * pc_y[k] * osf_388[k];

        t_583[k] = f_17 * nsf_309[k]
                   + f_3 * pc_y[k] * osf_389[k];

        t_584[k] = f_8 * nsf_299[k]
                   + f_1 * osd0_233[k]
                   - f_2 * osd1_233[k]
                   + f_3 * pc_z[k] * osf_389[k];
    }

#pragma omp simd aligned(t_585, t_586, t_587, pc_x, pc_y, pc_z, nsf_300, nsf_310, nsf_390, \
                         osd0_234, osd1_234, osf_390 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_585[k] = f_14 * nsf_390[k]
                   + f_1 * osd0_234[k]
                   - f_2 * osd1_234[k]
                   + f_3 * pc_x[k] * osf_390[k];

        t_586[k] = f_18 * nsf_310[k]
                   + f_3 * pc_y[k] * osf_390[k];

        t_587[k] = f_14 * nsf_300[k]
                   + f_3 * pc_z[k] * osf_390[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, pc_x, pc_y, nsf_312, nsf_393, nsf_395, osd0_237, \
                         osd0_239, osd1_237, osd1_239, osf_392, osf_393, \
                         osf_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = f_14 * nsf_393[k]
                   + f_4 * osd0_237[k]
                   - f_5 * osd1_237[k]
                   + f_3 * pc_x[k] * osf_393[k];

        t_589[k] = f_18 * nsf_312[k]
                   + f_3 * pc_y[k] * osf_392[k];

        t_590[k] = f_14 * nsf_395[k]
                   + f_4 * osd0_239[k]
                   - f_5 * osd1_239[k]
                   + f_3 * pc_x[k] * osf_395[k];
    }

#pragma omp simd aligned(t_591, t_592, t_593, t_594, pc_x, nsf_396, nsf_397, nsf_398, nsf_399, \
                         osf_396, osf_397, osf_398, osf_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_591[k] = f_14 * nsf_396[k]
                   + f_3 * pc_x[k] * osf_396[k];

        t_592[k] = f_14 * nsf_397[k]
                   + f_3 * pc_x[k] * osf_397[k];

        t_593[k] = f_14 * nsf_398[k]
                   + f_3 * pc_x[k] * osf_398[k];

        t_594[k] = f_14 * nsf_399[k]
                   + f_3 * pc_x[k] * osf_399[k];
    }

#pragma omp simd aligned(t_595, t_596, t_597, pc_y, pc_z, nsf_306, nsf_316, nsf_318, osd0_237, \
                         osd0_239, osd1_237, osd1_239, osf_396, \
                         osf_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_595[k] = f_18 * nsf_316[k]
                   + f_1 * osd0_237[k]
                   - f_2 * osd1_237[k]
                   + f_3 * pc_y[k] * osf_396[k];

        t_596[k] = f_14 * nsf_306[k]
                   + f_3 * pc_z[k] * osf_396[k];

        t_597[k] = f_18 * nsf_318[k]
                   + f_4 * osd0_239[k]
                   - f_5 * osd1_239[k]
                   + f_3 * pc_y[k] * osf_398[k];
    }

#pragma omp simd aligned(t_598, t_599, t_600, pc_x, pc_y, pc_z, nsf_309, nsf_319, nsf_400, \
                         osd0_239, osd0_240, osd1_239, osd1_240, osf_399, \
                         osf_400 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_598[k] = f_18 * nsf_319[k]
                   + f_3 * pc_y[k] * osf_399[k];

        t_599[k] = f_14 * nsf_309[k]
                   + f_1 * osd0_239[k]
                   - f_2 * osd1_239[k]
                   + f_3 * pc_z[k] * osf_399[k];

        t_600[k] = f_14 * nsf_400[k]
                   + f_1 * osd0_240[k]
                   - f_2 * osd1_240[k]
                   + f_3 * pc_x[k] * osf_400[k];
    }
}

static auto
compute_prim_osg_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsg0,
                                                          const size_t nsf, const size_t nsg1,
                                                          const size_t osd0, const size_t osd1,
                                                          const size_t osf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 4.5 / q;
    const auto f_13 = 4.0 / q;
    const auto f_14 = 1.5 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.5 / q;

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
    auto *t_718 = buffer.data(target + 718);
    auto *t_719 = buffer.data(target + 719);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsg0_525 = buffer.data(nsg0 + 525);
    const auto *nsg0_528 = buffer.data(nsg0 + 528);
    const auto *nsg0_530 = buffer.data(nsg0 + 530);
    const auto *nsg0_539 = buffer.data(nsg0 + 539);
    const auto *nsg0_540 = buffer.data(nsg0 + 540);
    const auto *nsg0_543 = buffer.data(nsg0 + 543);
    const auto *nsg0_550 = buffer.data(nsg0 + 550);

    const auto *nsf_310 = buffer.data(nsf + 310);
    const auto *nsf_316 = buffer.data(nsf + 316);
    const auto *nsf_319 = buffer.data(nsf + 319);
    const auto *nsf_320 = buffer.data(nsf + 320);
    const auto *nsf_322 = buffer.data(nsf + 322);
    const auto *nsf_326 = buffer.data(nsf + 326);
    const auto *nsf_328 = buffer.data(nsf + 328);
    const auto *nsf_329 = buffer.data(nsf + 329);
    const auto *nsf_330 = buffer.data(nsf + 330);
    const auto *nsf_332 = buffer.data(nsf + 332);
    const auto *nsf_336 = buffer.data(nsf + 336);
    const auto *nsf_338 = buffer.data(nsf + 338);
    const auto *nsf_339 = buffer.data(nsf + 339);
    const auto *nsf_340 = buffer.data(nsf + 340);
    const auto *nsf_342 = buffer.data(nsf + 342);
    const auto *nsf_346 = buffer.data(nsf + 346);
    const auto *nsf_348 = buffer.data(nsf + 348);
    const auto *nsf_349 = buffer.data(nsf + 349);
    const auto *nsf_350 = buffer.data(nsf + 350);
    const auto *nsf_351 = buffer.data(nsf + 351);
    const auto *nsf_352 = buffer.data(nsf + 352);
    const auto *nsf_356 = buffer.data(nsf + 356);
    const auto *nsf_358 = buffer.data(nsf + 358);
    const auto *nsf_359 = buffer.data(nsf + 359);
    const auto *nsf_360 = buffer.data(nsf + 360);
    const auto *nsf_366 = buffer.data(nsf + 366);
    const auto *nsf_369 = buffer.data(nsf + 369);
    const auto *nsf_370 = buffer.data(nsf + 370);
    const auto *nsf_372 = buffer.data(nsf + 372);
    const auto *nsf_376 = buffer.data(nsf + 376);
    const auto *nsf_378 = buffer.data(nsf + 378);
    const auto *nsf_379 = buffer.data(nsf + 379);
    const auto *nsf_380 = buffer.data(nsf + 380);
    const auto *nsf_382 = buffer.data(nsf + 382);
    const auto *nsf_386 = buffer.data(nsf + 386);
    const auto *nsf_388 = buffer.data(nsf + 388);
    const auto *nsf_389 = buffer.data(nsf + 389);
    const auto *nsf_403 = buffer.data(nsf + 403);
    const auto *nsf_405 = buffer.data(nsf + 405);
    const auto *nsf_406 = buffer.data(nsf + 406);
    const auto *nsf_407 = buffer.data(nsf + 407);
    const auto *nsf_408 = buffer.data(nsf + 408);
    const auto *nsf_409 = buffer.data(nsf + 409);
    const auto *nsf_410 = buffer.data(nsf + 410);
    const auto *nsf_413 = buffer.data(nsf + 413);
    const auto *nsf_415 = buffer.data(nsf + 415);
    const auto *nsf_416 = buffer.data(nsf + 416);
    const auto *nsf_417 = buffer.data(nsf + 417);
    const auto *nsf_418 = buffer.data(nsf + 418);
    const auto *nsf_419 = buffer.data(nsf + 419);
    const auto *nsf_420 = buffer.data(nsf + 420);
    const auto *nsf_423 = buffer.data(nsf + 423);
    const auto *nsf_425 = buffer.data(nsf + 425);
    const auto *nsf_426 = buffer.data(nsf + 426);
    const auto *nsf_427 = buffer.data(nsf + 427);
    const auto *nsf_428 = buffer.data(nsf + 428);
    const auto *nsf_429 = buffer.data(nsf + 429);
    const auto *nsf_436 = buffer.data(nsf + 436);
    const auto *nsf_437 = buffer.data(nsf + 437);
    const auto *nsf_438 = buffer.data(nsf + 438);
    const auto *nsf_439 = buffer.data(nsf + 439);
    const auto *nsf_440 = buffer.data(nsf + 440);
    const auto *nsf_445 = buffer.data(nsf + 445);
    const auto *nsf_446 = buffer.data(nsf + 446);
    const auto *nsf_447 = buffer.data(nsf + 447);
    const auto *nsf_449 = buffer.data(nsf + 449);
    const auto *nsf_450 = buffer.data(nsf + 450);
    const auto *nsf_453 = buffer.data(nsf + 453);
    const auto *nsf_456 = buffer.data(nsf + 456);
    const auto *nsf_458 = buffer.data(nsf + 458);
    const auto *nsf_459 = buffer.data(nsf + 459);
    const auto *nsf_465 = buffer.data(nsf + 465);
    const auto *nsf_466 = buffer.data(nsf + 466);
    const auto *nsf_467 = buffer.data(nsf + 467);
    const auto *nsf_468 = buffer.data(nsf + 468);
    const auto *nsf_469 = buffer.data(nsf + 469);
    const auto *nsf_470 = buffer.data(nsf + 470);
    const auto *nsf_473 = buffer.data(nsf + 473);
    const auto *nsf_475 = buffer.data(nsf + 475);
    const auto *nsf_476 = buffer.data(nsf + 476);
    const auto *nsf_477 = buffer.data(nsf + 477);
    const auto *nsf_478 = buffer.data(nsf + 478);
    const auto *nsf_479 = buffer.data(nsf + 479);

    const auto *nsg1_525 = buffer.data(nsg1 + 525);
    const auto *nsg1_528 = buffer.data(nsg1 + 528);
    const auto *nsg1_530 = buffer.data(nsg1 + 530);
    const auto *nsg1_539 = buffer.data(nsg1 + 539);
    const auto *nsg1_540 = buffer.data(nsg1 + 540);
    const auto *nsg1_543 = buffer.data(nsg1 + 543);
    const auto *nsg1_550 = buffer.data(nsg1 + 550);

    const auto *osd0_243 = buffer.data(osd0 + 243);
    const auto *osd0_245 = buffer.data(osd0 + 245);
    const auto *osd0_246 = buffer.data(osd0 + 246);
    const auto *osd0_249 = buffer.data(osd0 + 249);
    const auto *osd0_251 = buffer.data(osd0 + 251);
    const auto *osd0_252 = buffer.data(osd0 + 252);
    const auto *osd0_255 = buffer.data(osd0 + 255);
    const auto *osd0_257 = buffer.data(osd0 + 257);
    const auto *osd0_261 = buffer.data(osd0 + 261);
    const auto *osd0_263 = buffer.data(osd0 + 263);
    const auto *osd0_264 = buffer.data(osd0 + 264);
    const auto *osd0_267 = buffer.data(osd0 + 267);
    const auto *osd0_268 = buffer.data(osd0 + 268);
    const auto *osd0_269 = buffer.data(osd0 + 269);
    const auto *osd0_270 = buffer.data(osd0 + 270);
    const auto *osd0_273 = buffer.data(osd0 + 273);
    const auto *osd0_275 = buffer.data(osd0 + 275);
    const auto *osd0_281 = buffer.data(osd0 + 281);
    const auto *osd0_282 = buffer.data(osd0 + 282);
    const auto *osd0_285 = buffer.data(osd0 + 285);
    const auto *osd0_287 = buffer.data(osd0 + 287);

    const auto *osd1_243 = buffer.data(osd1 + 243);
    const auto *osd1_245 = buffer.data(osd1 + 245);
    const auto *osd1_246 = buffer.data(osd1 + 246);
    const auto *osd1_249 = buffer.data(osd1 + 249);
    const auto *osd1_251 = buffer.data(osd1 + 251);
    const auto *osd1_252 = buffer.data(osd1 + 252);
    const auto *osd1_255 = buffer.data(osd1 + 255);
    const auto *osd1_257 = buffer.data(osd1 + 257);
    const auto *osd1_261 = buffer.data(osd1 + 261);
    const auto *osd1_263 = buffer.data(osd1 + 263);
    const auto *osd1_264 = buffer.data(osd1 + 264);
    const auto *osd1_267 = buffer.data(osd1 + 267);
    const auto *osd1_268 = buffer.data(osd1 + 268);
    const auto *osd1_269 = buffer.data(osd1 + 269);
    const auto *osd1_270 = buffer.data(osd1 + 270);
    const auto *osd1_273 = buffer.data(osd1 + 273);
    const auto *osd1_275 = buffer.data(osd1 + 275);
    const auto *osd1_281 = buffer.data(osd1 + 281);
    const auto *osd1_282 = buffer.data(osd1 + 282);
    const auto *osd1_285 = buffer.data(osd1 + 285);
    const auto *osd1_287 = buffer.data(osd1 + 287);

    const auto *osf_400 = buffer.data(osf + 400);
    const auto *osf_402 = buffer.data(osf + 402);
    const auto *osf_403 = buffer.data(osf + 403);
    const auto *osf_405 = buffer.data(osf + 405);
    const auto *osf_406 = buffer.data(osf + 406);
    const auto *osf_407 = buffer.data(osf + 407);
    const auto *osf_408 = buffer.data(osf + 408);
    const auto *osf_409 = buffer.data(osf + 409);
    const auto *osf_410 = buffer.data(osf + 410);
    const auto *osf_412 = buffer.data(osf + 412);
    const auto *osf_413 = buffer.data(osf + 413);
    const auto *osf_415 = buffer.data(osf + 415);
    const auto *osf_416 = buffer.data(osf + 416);
    const auto *osf_417 = buffer.data(osf + 417);
    const auto *osf_418 = buffer.data(osf + 418);
    const auto *osf_419 = buffer.data(osf + 419);
    const auto *osf_420 = buffer.data(osf + 420);
    const auto *osf_422 = buffer.data(osf + 422);
    const auto *osf_423 = buffer.data(osf + 423);
    const auto *osf_425 = buffer.data(osf + 425);
    const auto *osf_426 = buffer.data(osf + 426);
    const auto *osf_427 = buffer.data(osf + 427);
    const auto *osf_428 = buffer.data(osf + 428);
    const auto *osf_429 = buffer.data(osf + 429);
    const auto *osf_430 = buffer.data(osf + 430);
    const auto *osf_432 = buffer.data(osf + 432);
    const auto *osf_436 = buffer.data(osf + 436);
    const auto *osf_437 = buffer.data(osf + 437);
    const auto *osf_438 = buffer.data(osf + 438);
    const auto *osf_439 = buffer.data(osf + 439);
    const auto *osf_440 = buffer.data(osf + 440);
    const auto *osf_441 = buffer.data(osf + 441);
    const auto *osf_442 = buffer.data(osf + 442);
    const auto *osf_445 = buffer.data(osf + 445);
    const auto *osf_446 = buffer.data(osf + 446);
    const auto *osf_447 = buffer.data(osf + 447);
    const auto *osf_448 = buffer.data(osf + 448);
    const auto *osf_449 = buffer.data(osf + 449);
    const auto *osf_450 = buffer.data(osf + 450);
    const auto *osf_451 = buffer.data(osf + 451);
    const auto *osf_452 = buffer.data(osf + 452);
    const auto *osf_453 = buffer.data(osf + 453);
    const auto *osf_456 = buffer.data(osf + 456);
    const auto *osf_457 = buffer.data(osf + 457);
    const auto *osf_458 = buffer.data(osf + 458);
    const auto *osf_459 = buffer.data(osf + 459);
    const auto *osf_460 = buffer.data(osf + 460);
    const auto *osf_462 = buffer.data(osf + 462);
    const auto *osf_465 = buffer.data(osf + 465);
    const auto *osf_466 = buffer.data(osf + 466);
    const auto *osf_467 = buffer.data(osf + 467);
    const auto *osf_468 = buffer.data(osf + 468);
    const auto *osf_469 = buffer.data(osf + 469);
    const auto *osf_470 = buffer.data(osf + 470);
    const auto *osf_472 = buffer.data(osf + 472);
    const auto *osf_473 = buffer.data(osf + 473);
    const auto *osf_475 = buffer.data(osf + 475);
    const auto *osf_476 = buffer.data(osf + 476);
    const auto *osf_477 = buffer.data(osf + 477);
    const auto *osf_478 = buffer.data(osf + 478);
    const auto *osf_479 = buffer.data(osf + 479);

#pragma omp simd aligned(t_601, t_602, t_603, t_604, pc_x, pc_y, pc_z, nsf_310, nsf_320, \
                         nsf_322, nsf_403, osd0_243, osd1_243, osf_400, osf_402, \
                         osf_403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_601[k] = f_16 * nsf_320[k]
                   + f_3 * pc_y[k] * osf_400[k];

        t_602[k] = f_16 * nsf_310[k]
                   + f_3 * pc_z[k] * osf_400[k];

        t_603[k] = f_14 * nsf_403[k]
                   + f_4 * osd0_243[k]
                   - f_5 * osd1_243[k]
                   + f_3 * pc_x[k] * osf_403[k];

        t_604[k] = f_16 * nsf_322[k]
                   + f_3 * pc_y[k] * osf_402[k];
    }

#pragma omp simd aligned(t_605, t_606, t_607, t_608, pc_x, nsf_405, nsf_406, nsf_407, nsf_408, \
                         osd0_245, osd1_245, osf_405, osf_406, osf_407, \
                         osf_408 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_605[k] = f_14 * nsf_405[k]
                   + f_4 * osd0_245[k]
                   - f_5 * osd1_245[k]
                   + f_3 * pc_x[k] * osf_405[k];

        t_606[k] = f_14 * nsf_406[k]
                   + f_3 * pc_x[k] * osf_406[k];

        t_607[k] = f_14 * nsf_407[k]
                   + f_3 * pc_x[k] * osf_407[k];

        t_608[k] = f_14 * nsf_408[k]
                   + f_3 * pc_x[k] * osf_408[k];
    }

#pragma omp simd aligned(t_609, t_610, t_611, pc_x, pc_y, pc_z, nsf_316, nsf_326, nsf_409, \
                         osd0_243, osd1_243, osf_406, osf_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_609[k] = f_14 * nsf_409[k]
                   + f_3 * pc_x[k] * osf_409[k];

        t_610[k] = f_16 * nsf_326[k]
                   + f_1 * osd0_243[k]
                   - f_2 * osd1_243[k]
                   + f_3 * pc_y[k] * osf_406[k];

        t_611[k] = f_16 * nsf_316[k]
                   + f_3 * pc_z[k] * osf_406[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, pc_y, pc_z, nsf_319, nsf_328, nsf_329, osd0_245, \
                         osd1_245, osf_408, osf_409 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = f_16 * nsf_328[k]
                   + f_4 * osd0_245[k]
                   - f_5 * osd1_245[k]
                   + f_3 * pc_y[k] * osf_408[k];

        t_613[k] = f_16 * nsf_329[k]
                   + f_3 * pc_y[k] * osf_409[k];

        t_614[k] = f_16 * nsf_319[k]
                   + f_1 * osd0_245[k]
                   - f_2 * osd1_245[k]
                   + f_3 * pc_z[k] * osf_409[k];
    }

#pragma omp simd aligned(t_615, t_616, t_617, pc_x, pc_y, pc_z, nsf_320, nsf_330, nsf_410, \
                         osd0_246, osd1_246, osf_410 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_615[k] = f_14 * nsf_410[k]
                   + f_1 * osd0_246[k]
                   - f_2 * osd1_246[k]
                   + f_3 * pc_x[k] * osf_410[k];

        t_616[k] = f_14 * nsf_330[k]
                   + f_3 * pc_y[k] * osf_410[k];

        t_617[k] = f_18 * nsf_320[k]
                   + f_3 * pc_z[k] * osf_410[k];
    }

#pragma omp simd aligned(t_618, t_619, t_620, pc_x, pc_y, nsf_332, nsf_413, nsf_415, osd0_249, \
                         osd0_251, osd1_249, osd1_251, osf_412, osf_413, \
                         osf_415 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_618[k] = f_14 * nsf_413[k]
                   + f_4 * osd0_249[k]
                   - f_5 * osd1_249[k]
                   + f_3 * pc_x[k] * osf_413[k];

        t_619[k] = f_14 * nsf_332[k]
                   + f_3 * pc_y[k] * osf_412[k];

        t_620[k] = f_14 * nsf_415[k]
                   + f_4 * osd0_251[k]
                   - f_5 * osd1_251[k]
                   + f_3 * pc_x[k] * osf_415[k];
    }

#pragma omp simd aligned(t_621, t_622, t_623, t_624, pc_x, nsf_416, nsf_417, nsf_418, nsf_419, \
                         osf_416, osf_417, osf_418, osf_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_621[k] = f_14 * nsf_416[k]
                   + f_3 * pc_x[k] * osf_416[k];

        t_622[k] = f_14 * nsf_417[k]
                   + f_3 * pc_x[k] * osf_417[k];

        t_623[k] = f_14 * nsf_418[k]
                   + f_3 * pc_x[k] * osf_418[k];

        t_624[k] = f_14 * nsf_419[k]
                   + f_3 * pc_x[k] * osf_419[k];
    }

#pragma omp simd aligned(t_625, t_626, t_627, pc_y, pc_z, nsf_326, nsf_336, nsf_338, osd0_249, \
                         osd0_251, osd1_249, osd1_251, osf_416, \
                         osf_418 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_625[k] = f_14 * nsf_336[k]
                   + f_1 * osd0_249[k]
                   - f_2 * osd1_249[k]
                   + f_3 * pc_y[k] * osf_416[k];

        t_626[k] = f_18 * nsf_326[k]
                   + f_3 * pc_z[k] * osf_416[k];

        t_627[k] = f_14 * nsf_338[k]
                   + f_4 * osd0_251[k]
                   - f_5 * osd1_251[k]
                   + f_3 * pc_y[k] * osf_418[k];
    }

#pragma omp simd aligned(t_628, t_629, t_630, pc_x, pc_y, pc_z, nsf_329, nsf_339, nsf_420, \
                         osd0_251, osd0_252, osd1_251, osd1_252, osf_419, \
                         osf_420 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_628[k] = f_14 * nsf_339[k]
                   + f_3 * pc_y[k] * osf_419[k];

        t_629[k] = f_18 * nsf_329[k]
                   + f_1 * osd0_251[k]
                   - f_2 * osd1_251[k]
                   + f_3 * pc_z[k] * osf_419[k];

        t_630[k] = f_14 * nsf_420[k]
                   + f_1 * osd0_252[k]
                   - f_2 * osd1_252[k]
                   + f_3 * pc_x[k] * osf_420[k];
    }

#pragma omp simd aligned(t_631, t_632, t_633, t_634, pc_x, pc_y, pc_z, nsf_330, nsf_340, \
                         nsf_342, nsf_423, osd0_255, osd1_255, osf_420, osf_422, \
                         osf_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_631[k] = f_8 * nsf_340[k]
                   + f_3 * pc_y[k] * osf_420[k];

        t_632[k] = f_17 * nsf_330[k]
                   + f_3 * pc_z[k] * osf_420[k];

        t_633[k] = f_14 * nsf_423[k]
                   + f_4 * osd0_255[k]
                   - f_5 * osd1_255[k]
                   + f_3 * pc_x[k] * osf_423[k];

        t_634[k] = f_8 * nsf_342[k]
                   + f_3 * pc_y[k] * osf_422[k];
    }

#pragma omp simd aligned(t_635, t_636, t_637, t_638, pc_x, nsf_425, nsf_426, nsf_427, nsf_428, \
                         osd0_257, osd1_257, osf_425, osf_426, osf_427, \
                         osf_428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_635[k] = f_14 * nsf_425[k]
                   + f_4 * osd0_257[k]
                   - f_5 * osd1_257[k]
                   + f_3 * pc_x[k] * osf_425[k];

        t_636[k] = f_14 * nsf_426[k]
                   + f_3 * pc_x[k] * osf_426[k];

        t_637[k] = f_14 * nsf_427[k]
                   + f_3 * pc_x[k] * osf_427[k];

        t_638[k] = f_14 * nsf_428[k]
                   + f_3 * pc_x[k] * osf_428[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, pc_x, pc_y, pc_z, nsf_336, nsf_346, nsf_429, \
                         osd0_255, osd1_255, osf_426, osf_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = f_14 * nsf_429[k]
                   + f_3 * pc_x[k] * osf_429[k];

        t_640[k] = f_8 * nsf_346[k]
                   + f_1 * osd0_255[k]
                   - f_2 * osd1_255[k]
                   + f_3 * pc_y[k] * osf_426[k];

        t_641[k] = f_17 * nsf_336[k]
                   + f_3 * pc_z[k] * osf_426[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, t_645, pa_y, pc_y, pc_z, nsg0_525, nsf_339, \
                         nsf_348, nsf_349, nsg1_525, osd0_257, osd1_257, osf_428, \
                         osf_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = f_8 * nsf_348[k]
                   + f_4 * osd0_257[k]
                   - f_5 * osd1_257[k]
                   + f_3 * pc_y[k] * osf_428[k];

        t_643[k] = f_8 * nsf_349[k]
                   + f_3 * pc_y[k] * osf_429[k];

        t_644[k] = f_17 * nsf_339[k]
                   + f_1 * osd0_257[k]
                   - f_2 * osd1_257[k]
                   + f_3 * pc_z[k] * osf_429[k];

        t_645[k] = pa_y[k] * nsg0_525[k]
                   - f_6 * pc_y[k] * nsg1_525[k];
    }

#pragma omp simd aligned(t_646, t_647, t_648, t_649, pa_y, pc_y, pc_z, nsg0_528, nsf_340, \
                         nsf_350, nsf_351, nsf_352, nsg1_528, osf_430, \
                         osf_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_646[k] = f_7 * nsf_350[k]
                   + f_3 * pc_y[k] * osf_430[k];

        t_647[k] = f_15 * nsf_340[k]
                   + f_3 * pc_z[k] * osf_430[k];

        t_648[k] = pa_y[k] * nsg0_528[k]
                   + f_8 * nsf_351[k]
                   - f_6 * pc_y[k] * nsg1_528[k];

        t_649[k] = f_7 * nsf_352[k]
                   + f_3 * pc_y[k] * osf_432[k];
    }

#pragma omp simd aligned(t_650, t_651, t_652, t_653, pa_y, pc_x, pc_y, nsg0_530, nsf_436, \
                         nsf_437, nsf_438, nsg1_530, osf_436, osf_437, \
                         osf_438 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_650[k] = pa_y[k] * nsg0_530[k]
                   - f_6 * pc_y[k] * nsg1_530[k];

        t_651[k] = f_14 * nsf_436[k]
                   + f_3 * pc_x[k] * osf_436[k];

        t_652[k] = f_14 * nsf_437[k]
                   + f_3 * pc_x[k] * osf_437[k];

        t_653[k] = f_14 * nsf_438[k]
                   + f_3 * pc_x[k] * osf_438[k];
    }

#pragma omp simd aligned(t_654, t_655, t_656, pc_x, pc_y, pc_z, nsf_346, nsf_356, nsf_439, \
                         osd0_261, osd1_261, osf_436, osf_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_654[k] = f_14 * nsf_439[k]
                   + f_3 * pc_x[k] * osf_439[k];

        t_655[k] = f_7 * nsf_356[k]
                   + f_1 * osd0_261[k]
                   - f_2 * osd1_261[k]
                   + f_3 * pc_y[k] * osf_436[k];

        t_656[k] = f_15 * nsf_346[k]
                   + f_3 * pc_z[k] * osf_436[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, pa_y, pc_y, nsg0_539, nsf_358, nsf_359, \
                         nsg1_539, osd0_263, osd1_263, osf_438, \
                         osf_439 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = f_7 * nsf_358[k]
                   + f_4 * osd0_263[k]
                   - f_5 * osd1_263[k]
                   + f_3 * pc_y[k] * osf_438[k];

        t_658[k] = f_7 * nsf_359[k]
                   + f_3 * pc_y[k] * osf_439[k];

        t_659[k] = pa_y[k] * nsg0_539[k]
                   - f_6 * pc_y[k] * nsg1_539[k];
    }

#pragma omp simd aligned(t_660, t_661, t_662, t_663, t_664, pc_x, pc_y, pc_z, nsf_350, \
                         nsf_440, osd0_264, osd1_264, osf_440, osf_441, \
                         osf_442 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_660[k] = f_14 * nsf_440[k]
                   + f_1 * osd0_264[k]
                   - f_2 * osd1_264[k]
                   + f_3 * pc_x[k] * osf_440[k];

        t_661[k] = f_3 * pc_y[k] * osf_440[k];

        t_662[k] = f_13 * nsf_350[k]
                   + f_3 * pc_z[k] * osf_440[k];

        t_663[k] = f_4 * osd0_264[k]
                   - f_5 * osd1_264[k]
                   + f_3 * pc_y[k] * osf_441[k];

        t_664[k] = f_3 * pc_y[k] * osf_442[k];
    }

#pragma omp simd aligned(t_665, t_666, t_667, t_668, pc_x, pc_y, nsf_445, nsf_446, nsf_447, \
                         osd0_269, osd1_269, osf_445, osf_446, \
                         osf_447 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_665[k] = f_14 * nsf_445[k]
                   + f_4 * osd0_269[k]
                   - f_5 * osd1_269[k]
                   + f_3 * pc_x[k] * osf_445[k];

        t_666[k] = f_14 * nsf_446[k]
                   + f_3 * pc_x[k] * osf_446[k];

        t_667[k] = f_14 * nsf_447[k]
                   + f_3 * pc_x[k] * osf_447[k];

        t_668[k] = f_3 * pc_y[k] * osf_445[k];
    }

#pragma omp simd aligned(t_669, t_670, t_671, pc_x, pc_y, nsf_449, osd0_267, osd0_268, \
                         osd1_267, osd1_268, osf_446, osf_447, \
                         osf_449 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_669[k] = f_14 * nsf_449[k]
                   + f_3 * pc_x[k] * osf_449[k];

        t_670[k] = f_1 * osd0_267[k]
                   - f_2 * osd1_267[k]
                   + f_3 * pc_y[k] * osf_446[k];

        t_671[k] = f_10 * osd0_268[k]
                   - f_11 * osd1_268[k]
                   + f_3 * pc_y[k] * osf_447[k];
    }

#pragma omp simd aligned(t_672, t_673, t_674, t_675, pc_x, pc_y, pc_z, nsf_359, nsf_450, \
                         osd0_269, osd0_270, osd1_269, osd1_270, osf_448, osf_449, \
                         osf_450 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_672[k] = f_4 * osd0_269[k]
                   - f_5 * osd1_269[k]
                   + f_3 * pc_y[k] * osf_448[k];

        t_673[k] = f_3 * pc_y[k] * osf_449[k];

        t_674[k] = f_13 * nsf_359[k]
                   + f_1 * osd0_269[k]
                   - f_2 * osd1_269[k]
                   + f_3 * pc_z[k] * osf_449[k];

        t_675[k] = f_8 * nsf_450[k]
                   + f_1 * osd0_270[k]
                   - f_2 * osd1_270[k]
                   + f_3 * pc_x[k] * osf_450[k];
    }

#pragma omp simd aligned(t_676, t_677, t_678, t_679, pc_x, pc_y, pc_z, nsf_360, nsf_453, \
                         osd0_273, osd1_273, osf_450, osf_451, \
                         osf_453 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_676[k] = f_12 * nsf_360[k]
                   + f_3 * pc_y[k] * osf_450[k];

        t_677[k] = f_3 * pc_z[k] * osf_450[k];

        t_678[k] = f_8 * nsf_453[k]
                   + f_4 * osd0_273[k]
                   - f_5 * osd1_273[k]
                   + f_3 * pc_x[k] * osf_453[k];

        t_679[k] = f_3 * pc_z[k] * osf_451[k];
    }

#pragma omp simd aligned(t_680, t_681, t_682, t_683, pc_x, pc_z, nsf_456, nsf_458, osd0_270, \
                         osd1_270, osf_452, osf_453, osf_456, osf_458 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_680[k] = f_4 * osd0_270[k]
                   - f_5 * osd1_270[k]
                   + f_3 * pc_z[k] * osf_452[k];

        t_681[k] = f_8 * nsf_456[k]
                   + f_3 * pc_x[k] * osf_456[k];

        t_682[k] = f_3 * pc_z[k] * osf_453[k];

        t_683[k] = f_8 * nsf_458[k]
                   + f_3 * pc_x[k] * osf_458[k];
    }

#pragma omp simd aligned(t_684, t_685, t_686, t_687, t_688, pc_x, pc_y, pc_z, nsf_366, \
                         nsf_369, nsf_459, osd0_273, osd1_273, osf_456, osf_457, \
                         osf_459 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_684[k] = f_8 * nsf_459[k]
                   + f_3 * pc_x[k] * osf_459[k];

        t_685[k] = f_12 * nsf_366[k]
                   + f_1 * osd0_273[k]
                   - f_2 * osd1_273[k]
                   + f_3 * pc_y[k] * osf_456[k];

        t_686[k] = f_3 * pc_z[k] * osf_456[k];

        t_687[k] = f_4 * osd0_273[k]
                   - f_5 * osd1_273[k]
                   + f_3 * pc_z[k] * osf_457[k];

        t_688[k] = f_12 * nsf_369[k]
                   + f_3 * pc_y[k] * osf_459[k];
    }

#pragma omp simd aligned(t_689, t_690, t_691, t_692, pa_z, pc_y, pc_z, nsg0_540, nsf_360, \
                         nsf_370, nsg1_540, osd0_275, osd1_275, osf_459, \
                         osf_460 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_689[k] = f_1 * osd0_275[k]
                   - f_2 * osd1_275[k]
                   + f_3 * pc_z[k] * osf_459[k];

        t_690[k] = pa_z[k] * nsg0_540[k]
                   - f_6 * pc_z[k] * nsg1_540[k];

        t_691[k] = f_13 * nsf_370[k]
                   + f_3 * pc_y[k] * osf_460[k];

        t_692[k] = f_7 * nsf_360[k]
                   + f_3 * pc_z[k] * osf_460[k];
    }

#pragma omp simd aligned(t_693, t_694, t_695, pa_z, pc_x, pc_y, pc_z, nsg0_543, nsf_372, \
                         nsf_465, nsg1_543, osd0_281, osd1_281, osf_462, \
                         osf_465 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_693[k] = pa_z[k] * nsg0_543[k]
                   - f_6 * pc_z[k] * nsg1_543[k];

        t_694[k] = f_13 * nsf_372[k]
                   + f_3 * pc_y[k] * osf_462[k];

        t_695[k] = f_8 * nsf_465[k]
                   + f_4 * osd0_281[k]
                   - f_5 * osd1_281[k]
                   + f_3 * pc_x[k] * osf_465[k];
    }

#pragma omp simd aligned(t_696, t_697, t_698, t_699, pc_x, nsf_466, nsf_467, nsf_468, nsf_469, \
                         osf_466, osf_467, osf_468, osf_469 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_696[k] = f_8 * nsf_466[k]
                   + f_3 * pc_x[k] * osf_466[k];

        t_697[k] = f_8 * nsf_467[k]
                   + f_3 * pc_x[k] * osf_467[k];

        t_698[k] = f_8 * nsf_468[k]
                   + f_3 * pc_x[k] * osf_468[k];

        t_699[k] = f_8 * nsf_469[k]
                   + f_3 * pc_x[k] * osf_469[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, pa_z, pc_y, pc_z, nsg0_550, nsf_366, nsf_378, \
                         nsg1_550, osd0_281, osd1_281, osf_466, \
                         osf_468 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = pa_z[k] * nsg0_550[k]
                   - f_6 * pc_z[k] * nsg1_550[k];

        t_701[k] = f_7 * nsf_366[k]
                   + f_3 * pc_z[k] * osf_466[k];

        t_702[k] = f_13 * nsf_378[k]
                   + f_4 * osd0_281[k]
                   - f_5 * osd1_281[k]
                   + f_3 * pc_y[k] * osf_468[k];
    }

#pragma omp simd aligned(t_703, t_704, t_705, pc_x, pc_y, pc_z, nsf_369, nsf_379, nsf_470, \
                         osd0_281, osd0_282, osd1_281, osd1_282, osf_469, \
                         osf_470 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_703[k] = f_13 * nsf_379[k]
                   + f_3 * pc_y[k] * osf_469[k];

        t_704[k] = f_7 * nsf_369[k]
                   + f_1 * osd0_281[k]
                   - f_2 * osd1_281[k]
                   + f_3 * pc_z[k] * osf_469[k];

        t_705[k] = f_8 * nsf_470[k]
                   + f_1 * osd0_282[k]
                   - f_2 * osd1_282[k]
                   + f_3 * pc_x[k] * osf_470[k];
    }

#pragma omp simd aligned(t_706, t_707, t_708, t_709, pc_x, pc_y, pc_z, nsf_370, nsf_380, \
                         nsf_382, nsf_473, osd0_285, osd1_285, osf_470, osf_472, \
                         osf_473 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_706[k] = f_15 * nsf_380[k]
                   + f_3 * pc_y[k] * osf_470[k];

        t_707[k] = f_8 * nsf_370[k]
                   + f_3 * pc_z[k] * osf_470[k];

        t_708[k] = f_8 * nsf_473[k]
                   + f_4 * osd0_285[k]
                   - f_5 * osd1_285[k]
                   + f_3 * pc_x[k] * osf_473[k];

        t_709[k] = f_15 * nsf_382[k]
                   + f_3 * pc_y[k] * osf_472[k];
    }

#pragma omp simd aligned(t_710, t_711, t_712, t_713, pc_x, nsf_475, nsf_476, nsf_477, nsf_478, \
                         osd0_287, osd1_287, osf_475, osf_476, osf_477, \
                         osf_478 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_710[k] = f_8 * nsf_475[k]
                   + f_4 * osd0_287[k]
                   - f_5 * osd1_287[k]
                   + f_3 * pc_x[k] * osf_475[k];

        t_711[k] = f_8 * nsf_476[k]
                   + f_3 * pc_x[k] * osf_476[k];

        t_712[k] = f_8 * nsf_477[k]
                   + f_3 * pc_x[k] * osf_477[k];

        t_713[k] = f_8 * nsf_478[k]
                   + f_3 * pc_x[k] * osf_478[k];
    }

#pragma omp simd aligned(t_714, t_715, t_716, pc_x, pc_y, pc_z, nsf_376, nsf_386, nsf_479, \
                         osd0_285, osd1_285, osf_476, osf_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_714[k] = f_8 * nsf_479[k]
                   + f_3 * pc_x[k] * osf_479[k];

        t_715[k] = f_15 * nsf_386[k]
                   + f_1 * osd0_285[k]
                   - f_2 * osd1_285[k]
                   + f_3 * pc_y[k] * osf_476[k];

        t_716[k] = f_8 * nsf_376[k]
                   + f_3 * pc_z[k] * osf_476[k];
    }

#pragma omp simd aligned(t_717, t_718, t_719, pc_y, pc_z, nsf_379, nsf_388, nsf_389, osd0_287, \
                         osd1_287, osf_478, osf_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_717[k] = f_15 * nsf_388[k]
                   + f_4 * osd0_287[k]
                   - f_5 * osd1_287[k]
                   + f_3 * pc_y[k] * osf_478[k];

        t_718[k] = f_15 * nsf_389[k]
                   + f_3 * pc_y[k] * osf_479[k];

        t_719[k] = f_8 * nsf_379[k]
                   + f_1 * osd0_287[k]
                   - f_2 * osd1_287[k]
                   + f_3 * pc_z[k] * osf_479[k];
    }
}

static auto
compute_prim_osg_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsg0,
                                                          const size_t nsf, const size_t nsg1,
                                                          const size_t osd0, const size_t osd1,
                                                          const size_t osf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 5.0 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 4.5 / q;
    const auto f_13 = 4.0 / q;
    const auto f_14 = 1.5 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.5 / q;

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
    auto *t_834 = buffer.data(target + 834);
    auto *t_835 = buffer.data(target + 835);
    auto *t_836 = buffer.data(target + 836);
    auto *t_837 = buffer.data(target + 837);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsg0_660 = buffer.data(nsg0 + 660);
    const auto *nsg0_663 = buffer.data(nsg0 + 663);
    const auto *nsg0_665 = buffer.data(nsg0 + 665);
    const auto *nsg0_674 = buffer.data(nsg0 + 674);
    const auto *nsg0_825 = buffer.data(nsg0 + 825);
    const auto *nsg0_828 = buffer.data(nsg0 + 828);
    const auto *nsg0_835 = buffer.data(nsg0 + 835);
    const auto *nsg0_837 = buffer.data(nsg0 + 837);

    const auto *nsf_380 = buffer.data(nsf + 380);
    const auto *nsf_386 = buffer.data(nsf + 386);
    const auto *nsf_389 = buffer.data(nsf + 389);
    const auto *nsf_390 = buffer.data(nsf + 390);
    const auto *nsf_392 = buffer.data(nsf + 392);
    const auto *nsf_396 = buffer.data(nsf + 396);
    const auto *nsf_398 = buffer.data(nsf + 398);
    const auto *nsf_399 = buffer.data(nsf + 399);
    const auto *nsf_400 = buffer.data(nsf + 400);
    const auto *nsf_402 = buffer.data(nsf + 402);
    const auto *nsf_406 = buffer.data(nsf + 406);
    const auto *nsf_408 = buffer.data(nsf + 408);
    const auto *nsf_409 = buffer.data(nsf + 409);
    const auto *nsf_410 = buffer.data(nsf + 410);
    const auto *nsf_412 = buffer.data(nsf + 412);
    const auto *nsf_416 = buffer.data(nsf + 416);
    const auto *nsf_418 = buffer.data(nsf + 418);
    const auto *nsf_419 = buffer.data(nsf + 419);
    const auto *nsf_420 = buffer.data(nsf + 420);
    const auto *nsf_422 = buffer.data(nsf + 422);
    const auto *nsf_426 = buffer.data(nsf + 426);
    const auto *nsf_428 = buffer.data(nsf + 428);
    const auto *nsf_429 = buffer.data(nsf + 429);
    const auto *nsf_430 = buffer.data(nsf + 430);
    const auto *nsf_432 = buffer.data(nsf + 432);
    const auto *nsf_436 = buffer.data(nsf + 436);
    const auto *nsf_438 = buffer.data(nsf + 438);
    const auto *nsf_439 = buffer.data(nsf + 439);
    const auto *nsf_440 = buffer.data(nsf + 440);
    const auto *nsf_441 = buffer.data(nsf + 441);
    const auto *nsf_442 = buffer.data(nsf + 442);
    const auto *nsf_446 = buffer.data(nsf + 446);
    const auto *nsf_448 = buffer.data(nsf + 448);
    const auto *nsf_449 = buffer.data(nsf + 449);
    const auto *nsf_450 = buffer.data(nsf + 450);
    const auto *nsf_480 = buffer.data(nsf + 480);
    const auto *nsf_483 = buffer.data(nsf + 483);
    const auto *nsf_485 = buffer.data(nsf + 485);
    const auto *nsf_486 = buffer.data(nsf + 486);
    const auto *nsf_487 = buffer.data(nsf + 487);
    const auto *nsf_488 = buffer.data(nsf + 488);
    const auto *nsf_489 = buffer.data(nsf + 489);
    const auto *nsf_490 = buffer.data(nsf + 490);
    const auto *nsf_493 = buffer.data(nsf + 493);
    const auto *nsf_495 = buffer.data(nsf + 495);
    const auto *nsf_496 = buffer.data(nsf + 496);
    const auto *nsf_497 = buffer.data(nsf + 497);
    const auto *nsf_498 = buffer.data(nsf + 498);
    const auto *nsf_499 = buffer.data(nsf + 499);
    const auto *nsf_500 = buffer.data(nsf + 500);
    const auto *nsf_503 = buffer.data(nsf + 503);
    const auto *nsf_505 = buffer.data(nsf + 505);
    const auto *nsf_506 = buffer.data(nsf + 506);
    const auto *nsf_507 = buffer.data(nsf + 507);
    const auto *nsf_508 = buffer.data(nsf + 508);
    const auto *nsf_509 = buffer.data(nsf + 509);
    const auto *nsf_510 = buffer.data(nsf + 510);
    const auto *nsf_513 = buffer.data(nsf + 513);
    const auto *nsf_515 = buffer.data(nsf + 515);
    const auto *nsf_516 = buffer.data(nsf + 516);
    const auto *nsf_517 = buffer.data(nsf + 517);
    const auto *nsf_518 = buffer.data(nsf + 518);
    const auto *nsf_519 = buffer.data(nsf + 519);
    const auto *nsf_520 = buffer.data(nsf + 520);
    const auto *nsf_523 = buffer.data(nsf + 523);
    const auto *nsf_525 = buffer.data(nsf + 525);
    const auto *nsf_526 = buffer.data(nsf + 526);
    const auto *nsf_527 = buffer.data(nsf + 527);
    const auto *nsf_528 = buffer.data(nsf + 528);
    const auto *nsf_529 = buffer.data(nsf + 529);
    const auto *nsf_536 = buffer.data(nsf + 536);
    const auto *nsf_537 = buffer.data(nsf + 537);
    const auto *nsf_538 = buffer.data(nsf + 538);
    const auto *nsf_539 = buffer.data(nsf + 539);
    const auto *nsf_540 = buffer.data(nsf + 540);
    const auto *nsf_545 = buffer.data(nsf + 545);
    const auto *nsf_546 = buffer.data(nsf + 546);
    const auto *nsf_547 = buffer.data(nsf + 547);
    const auto *nsf_549 = buffer.data(nsf + 549);
    const auto *nsf_550 = buffer.data(nsf + 550);
    const auto *nsf_553 = buffer.data(nsf + 553);
    const auto *nsf_556 = buffer.data(nsf + 556);
    const auto *nsf_558 = buffer.data(nsf + 558);
    const auto *nsf_559 = buffer.data(nsf + 559);

    const auto *nsg1_660 = buffer.data(nsg1 + 660);
    const auto *nsg1_663 = buffer.data(nsg1 + 663);
    const auto *nsg1_665 = buffer.data(nsg1 + 665);
    const auto *nsg1_674 = buffer.data(nsg1 + 674);
    const auto *nsg1_825 = buffer.data(nsg1 + 825);
    const auto *nsg1_828 = buffer.data(nsg1 + 828);
    const auto *nsg1_835 = buffer.data(nsg1 + 835);
    const auto *nsg1_837 = buffer.data(nsg1 + 837);

    const auto *osd0_288 = buffer.data(osd0 + 288);
    const auto *osd0_291 = buffer.data(osd0 + 291);
    const auto *osd0_293 = buffer.data(osd0 + 293);
    const auto *osd0_294 = buffer.data(osd0 + 294);
    const auto *osd0_297 = buffer.data(osd0 + 297);
    const auto *osd0_299 = buffer.data(osd0 + 299);
    const auto *osd0_300 = buffer.data(osd0 + 300);
    const auto *osd0_303 = buffer.data(osd0 + 303);
    const auto *osd0_305 = buffer.data(osd0 + 305);
    const auto *osd0_306 = buffer.data(osd0 + 306);
    const auto *osd0_309 = buffer.data(osd0 + 309);
    const auto *osd0_311 = buffer.data(osd0 + 311);
    const auto *osd0_312 = buffer.data(osd0 + 312);
    const auto *osd0_315 = buffer.data(osd0 + 315);
    const auto *osd0_317 = buffer.data(osd0 + 317);
    const auto *osd0_321 = buffer.data(osd0 + 321);
    const auto *osd0_323 = buffer.data(osd0 + 323);
    const auto *osd0_324 = buffer.data(osd0 + 324);
    const auto *osd0_327 = buffer.data(osd0 + 327);
    const auto *osd0_328 = buffer.data(osd0 + 328);
    const auto *osd0_329 = buffer.data(osd0 + 329);
    const auto *osd0_330 = buffer.data(osd0 + 330);

    const auto *osd1_288 = buffer.data(osd1 + 288);
    const auto *osd1_291 = buffer.data(osd1 + 291);
    const auto *osd1_293 = buffer.data(osd1 + 293);
    const auto *osd1_294 = buffer.data(osd1 + 294);
    const auto *osd1_297 = buffer.data(osd1 + 297);
    const auto *osd1_299 = buffer.data(osd1 + 299);
    const auto *osd1_300 = buffer.data(osd1 + 300);
    const auto *osd1_303 = buffer.data(osd1 + 303);
    const auto *osd1_305 = buffer.data(osd1 + 305);
    const auto *osd1_306 = buffer.data(osd1 + 306);
    const auto *osd1_309 = buffer.data(osd1 + 309);
    const auto *osd1_311 = buffer.data(osd1 + 311);
    const auto *osd1_312 = buffer.data(osd1 + 312);
    const auto *osd1_315 = buffer.data(osd1 + 315);
    const auto *osd1_317 = buffer.data(osd1 + 317);
    const auto *osd1_321 = buffer.data(osd1 + 321);
    const auto *osd1_323 = buffer.data(osd1 + 323);
    const auto *osd1_324 = buffer.data(osd1 + 324);
    const auto *osd1_327 = buffer.data(osd1 + 327);
    const auto *osd1_328 = buffer.data(osd1 + 328);
    const auto *osd1_329 = buffer.data(osd1 + 329);
    const auto *osd1_330 = buffer.data(osd1 + 330);

    const auto *osf_480 = buffer.data(osf + 480);
    const auto *osf_482 = buffer.data(osf + 482);
    const auto *osf_483 = buffer.data(osf + 483);
    const auto *osf_485 = buffer.data(osf + 485);
    const auto *osf_486 = buffer.data(osf + 486);
    const auto *osf_487 = buffer.data(osf + 487);
    const auto *osf_488 = buffer.data(osf + 488);
    const auto *osf_489 = buffer.data(osf + 489);
    const auto *osf_490 = buffer.data(osf + 490);
    const auto *osf_492 = buffer.data(osf + 492);
    const auto *osf_493 = buffer.data(osf + 493);
    const auto *osf_495 = buffer.data(osf + 495);
    const auto *osf_496 = buffer.data(osf + 496);
    const auto *osf_497 = buffer.data(osf + 497);
    const auto *osf_498 = buffer.data(osf + 498);
    const auto *osf_499 = buffer.data(osf + 499);
    const auto *osf_500 = buffer.data(osf + 500);
    const auto *osf_502 = buffer.data(osf + 502);
    const auto *osf_503 = buffer.data(osf + 503);
    const auto *osf_505 = buffer.data(osf + 505);
    const auto *osf_506 = buffer.data(osf + 506);
    const auto *osf_507 = buffer.data(osf + 507);
    const auto *osf_508 = buffer.data(osf + 508);
    const auto *osf_509 = buffer.data(osf + 509);
    const auto *osf_510 = buffer.data(osf + 510);
    const auto *osf_512 = buffer.data(osf + 512);
    const auto *osf_513 = buffer.data(osf + 513);
    const auto *osf_515 = buffer.data(osf + 515);
    const auto *osf_516 = buffer.data(osf + 516);
    const auto *osf_517 = buffer.data(osf + 517);
    const auto *osf_518 = buffer.data(osf + 518);
    const auto *osf_519 = buffer.data(osf + 519);
    const auto *osf_520 = buffer.data(osf + 520);
    const auto *osf_522 = buffer.data(osf + 522);
    const auto *osf_523 = buffer.data(osf + 523);
    const auto *osf_525 = buffer.data(osf + 525);
    const auto *osf_526 = buffer.data(osf + 526);
    const auto *osf_527 = buffer.data(osf + 527);
    const auto *osf_528 = buffer.data(osf + 528);
    const auto *osf_529 = buffer.data(osf + 529);
    const auto *osf_530 = buffer.data(osf + 530);
    const auto *osf_532 = buffer.data(osf + 532);
    const auto *osf_536 = buffer.data(osf + 536);
    const auto *osf_537 = buffer.data(osf + 537);
    const auto *osf_538 = buffer.data(osf + 538);
    const auto *osf_539 = buffer.data(osf + 539);
    const auto *osf_540 = buffer.data(osf + 540);
    const auto *osf_541 = buffer.data(osf + 541);
    const auto *osf_542 = buffer.data(osf + 542);
    const auto *osf_545 = buffer.data(osf + 545);
    const auto *osf_546 = buffer.data(osf + 546);
    const auto *osf_547 = buffer.data(osf + 547);
    const auto *osf_548 = buffer.data(osf + 548);
    const auto *osf_549 = buffer.data(osf + 549);
    const auto *osf_550 = buffer.data(osf + 550);
    const auto *osf_551 = buffer.data(osf + 551);
    const auto *osf_552 = buffer.data(osf + 552);
    const auto *osf_553 = buffer.data(osf + 553);
    const auto *osf_556 = buffer.data(osf + 556);
    const auto *osf_558 = buffer.data(osf + 558);
    const auto *osf_559 = buffer.data(osf + 559);

#pragma omp simd aligned(t_720, t_721, t_722, pc_x, pc_y, pc_z, nsf_380, nsf_390, nsf_480, \
                         osd0_288, osd1_288, osf_480 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_720[k] = f_8 * nsf_480[k]
                   + f_1 * osd0_288[k]
                   - f_2 * osd1_288[k]
                   + f_3 * pc_x[k] * osf_480[k];

        t_721[k] = f_17 * nsf_390[k]
                   + f_3 * pc_y[k] * osf_480[k];

        t_722[k] = f_14 * nsf_380[k]
                   + f_3 * pc_z[k] * osf_480[k];
    }

#pragma omp simd aligned(t_723, t_724, t_725, pc_x, pc_y, nsf_392, nsf_483, nsf_485, osd0_291, \
                         osd0_293, osd1_291, osd1_293, osf_482, osf_483, \
                         osf_485 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_723[k] = f_8 * nsf_483[k]
                   + f_4 * osd0_291[k]
                   - f_5 * osd1_291[k]
                   + f_3 * pc_x[k] * osf_483[k];

        t_724[k] = f_17 * nsf_392[k]
                   + f_3 * pc_y[k] * osf_482[k];

        t_725[k] = f_8 * nsf_485[k]
                   + f_4 * osd0_293[k]
                   - f_5 * osd1_293[k]
                   + f_3 * pc_x[k] * osf_485[k];
    }

#pragma omp simd aligned(t_726, t_727, t_728, t_729, pc_x, nsf_486, nsf_487, nsf_488, nsf_489, \
                         osf_486, osf_487, osf_488, osf_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_726[k] = f_8 * nsf_486[k]
                   + f_3 * pc_x[k] * osf_486[k];

        t_727[k] = f_8 * nsf_487[k]
                   + f_3 * pc_x[k] * osf_487[k];

        t_728[k] = f_8 * nsf_488[k]
                   + f_3 * pc_x[k] * osf_488[k];

        t_729[k] = f_8 * nsf_489[k]
                   + f_3 * pc_x[k] * osf_489[k];
    }

#pragma omp simd aligned(t_730, t_731, t_732, pc_y, pc_z, nsf_386, nsf_396, nsf_398, osd0_291, \
                         osd0_293, osd1_291, osd1_293, osf_486, \
                         osf_488 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_730[k] = f_17 * nsf_396[k]
                   + f_1 * osd0_291[k]
                   - f_2 * osd1_291[k]
                   + f_3 * pc_y[k] * osf_486[k];

        t_731[k] = f_14 * nsf_386[k]
                   + f_3 * pc_z[k] * osf_486[k];

        t_732[k] = f_17 * nsf_398[k]
                   + f_4 * osd0_293[k]
                   - f_5 * osd1_293[k]
                   + f_3 * pc_y[k] * osf_488[k];
    }

#pragma omp simd aligned(t_733, t_734, t_735, pc_x, pc_y, pc_z, nsf_389, nsf_399, nsf_490, \
                         osd0_293, osd0_294, osd1_293, osd1_294, osf_489, \
                         osf_490 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_733[k] = f_17 * nsf_399[k]
                   + f_3 * pc_y[k] * osf_489[k];

        t_734[k] = f_14 * nsf_389[k]
                   + f_1 * osd0_293[k]
                   - f_2 * osd1_293[k]
                   + f_3 * pc_z[k] * osf_489[k];

        t_735[k] = f_8 * nsf_490[k]
                   + f_1 * osd0_294[k]
                   - f_2 * osd1_294[k]
                   + f_3 * pc_x[k] * osf_490[k];
    }

#pragma omp simd aligned(t_736, t_737, t_738, t_739, pc_x, pc_y, pc_z, nsf_390, nsf_400, \
                         nsf_402, nsf_493, osd0_297, osd1_297, osf_490, osf_492, \
                         osf_493 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_736[k] = f_18 * nsf_400[k]
                   + f_3 * pc_y[k] * osf_490[k];

        t_737[k] = f_16 * nsf_390[k]
                   + f_3 * pc_z[k] * osf_490[k];

        t_738[k] = f_8 * nsf_493[k]
                   + f_4 * osd0_297[k]
                   - f_5 * osd1_297[k]
                   + f_3 * pc_x[k] * osf_493[k];

        t_739[k] = f_18 * nsf_402[k]
                   + f_3 * pc_y[k] * osf_492[k];
    }

#pragma omp simd aligned(t_740, t_741, t_742, t_743, pc_x, nsf_495, nsf_496, nsf_497, nsf_498, \
                         osd0_299, osd1_299, osf_495, osf_496, osf_497, \
                         osf_498 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_740[k] = f_8 * nsf_495[k]
                   + f_4 * osd0_299[k]
                   - f_5 * osd1_299[k]
                   + f_3 * pc_x[k] * osf_495[k];

        t_741[k] = f_8 * nsf_496[k]
                   + f_3 * pc_x[k] * osf_496[k];

        t_742[k] = f_8 * nsf_497[k]
                   + f_3 * pc_x[k] * osf_497[k];

        t_743[k] = f_8 * nsf_498[k]
                   + f_3 * pc_x[k] * osf_498[k];
    }

#pragma omp simd aligned(t_744, t_745, t_746, pc_x, pc_y, pc_z, nsf_396, nsf_406, nsf_499, \
                         osd0_297, osd1_297, osf_496, osf_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_744[k] = f_8 * nsf_499[k]
                   + f_3 * pc_x[k] * osf_499[k];

        t_745[k] = f_18 * nsf_406[k]
                   + f_1 * osd0_297[k]
                   - f_2 * osd1_297[k]
                   + f_3 * pc_y[k] * osf_496[k];

        t_746[k] = f_16 * nsf_396[k]
                   + f_3 * pc_z[k] * osf_496[k];
    }

#pragma omp simd aligned(t_747, t_748, t_749, pc_y, pc_z, nsf_399, nsf_408, nsf_409, osd0_299, \
                         osd1_299, osf_498, osf_499 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_747[k] = f_18 * nsf_408[k]
                   + f_4 * osd0_299[k]
                   - f_5 * osd1_299[k]
                   + f_3 * pc_y[k] * osf_498[k];

        t_748[k] = f_18 * nsf_409[k]
                   + f_3 * pc_y[k] * osf_499[k];

        t_749[k] = f_16 * nsf_399[k]
                   + f_1 * osd0_299[k]
                   - f_2 * osd1_299[k]
                   + f_3 * pc_z[k] * osf_499[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, pc_x, pc_y, pc_z, nsf_400, nsf_410, nsf_500, \
                         osd0_300, osd1_300, osf_500 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = f_8 * nsf_500[k]
                   + f_1 * osd0_300[k]
                   - f_2 * osd1_300[k]
                   + f_3 * pc_x[k] * osf_500[k];

        t_751[k] = f_16 * nsf_410[k]
                   + f_3 * pc_y[k] * osf_500[k];

        t_752[k] = f_18 * nsf_400[k]
                   + f_3 * pc_z[k] * osf_500[k];
    }

#pragma omp simd aligned(t_753, t_754, t_755, pc_x, pc_y, nsf_412, nsf_503, nsf_505, osd0_303, \
                         osd0_305, osd1_303, osd1_305, osf_502, osf_503, \
                         osf_505 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_753[k] = f_8 * nsf_503[k]
                   + f_4 * osd0_303[k]
                   - f_5 * osd1_303[k]
                   + f_3 * pc_x[k] * osf_503[k];

        t_754[k] = f_16 * nsf_412[k]
                   + f_3 * pc_y[k] * osf_502[k];

        t_755[k] = f_8 * nsf_505[k]
                   + f_4 * osd0_305[k]
                   - f_5 * osd1_305[k]
                   + f_3 * pc_x[k] * osf_505[k];
    }

#pragma omp simd aligned(t_756, t_757, t_758, t_759, pc_x, nsf_506, nsf_507, nsf_508, nsf_509, \
                         osf_506, osf_507, osf_508, osf_509 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_756[k] = f_8 * nsf_506[k]
                   + f_3 * pc_x[k] * osf_506[k];

        t_757[k] = f_8 * nsf_507[k]
                   + f_3 * pc_x[k] * osf_507[k];

        t_758[k] = f_8 * nsf_508[k]
                   + f_3 * pc_x[k] * osf_508[k];

        t_759[k] = f_8 * nsf_509[k]
                   + f_3 * pc_x[k] * osf_509[k];
    }

#pragma omp simd aligned(t_760, t_761, t_762, pc_y, pc_z, nsf_406, nsf_416, nsf_418, osd0_303, \
                         osd0_305, osd1_303, osd1_305, osf_506, \
                         osf_508 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_760[k] = f_16 * nsf_416[k]
                   + f_1 * osd0_303[k]
                   - f_2 * osd1_303[k]
                   + f_3 * pc_y[k] * osf_506[k];

        t_761[k] = f_18 * nsf_406[k]
                   + f_3 * pc_z[k] * osf_506[k];

        t_762[k] = f_16 * nsf_418[k]
                   + f_4 * osd0_305[k]
                   - f_5 * osd1_305[k]
                   + f_3 * pc_y[k] * osf_508[k];
    }

#pragma omp simd aligned(t_763, t_764, t_765, pc_x, pc_y, pc_z, nsf_409, nsf_419, nsf_510, \
                         osd0_305, osd0_306, osd1_305, osd1_306, osf_509, \
                         osf_510 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_763[k] = f_16 * nsf_419[k]
                   + f_3 * pc_y[k] * osf_509[k];

        t_764[k] = f_18 * nsf_409[k]
                   + f_1 * osd0_305[k]
                   - f_2 * osd1_305[k]
                   + f_3 * pc_z[k] * osf_509[k];

        t_765[k] = f_8 * nsf_510[k]
                   + f_1 * osd0_306[k]
                   - f_2 * osd1_306[k]
                   + f_3 * pc_x[k] * osf_510[k];
    }

#pragma omp simd aligned(t_766, t_767, t_768, t_769, pc_x, pc_y, pc_z, nsf_410, nsf_420, \
                         nsf_422, nsf_513, osd0_309, osd1_309, osf_510, osf_512, \
                         osf_513 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_766[k] = f_14 * nsf_420[k]
                   + f_3 * pc_y[k] * osf_510[k];

        t_767[k] = f_17 * nsf_410[k]
                   + f_3 * pc_z[k] * osf_510[k];

        t_768[k] = f_8 * nsf_513[k]
                   + f_4 * osd0_309[k]
                   - f_5 * osd1_309[k]
                   + f_3 * pc_x[k] * osf_513[k];

        t_769[k] = f_14 * nsf_422[k]
                   + f_3 * pc_y[k] * osf_512[k];
    }

#pragma omp simd aligned(t_770, t_771, t_772, t_773, pc_x, nsf_515, nsf_516, nsf_517, nsf_518, \
                         osd0_311, osd1_311, osf_515, osf_516, osf_517, \
                         osf_518 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_770[k] = f_8 * nsf_515[k]
                   + f_4 * osd0_311[k]
                   - f_5 * osd1_311[k]
                   + f_3 * pc_x[k] * osf_515[k];

        t_771[k] = f_8 * nsf_516[k]
                   + f_3 * pc_x[k] * osf_516[k];

        t_772[k] = f_8 * nsf_517[k]
                   + f_3 * pc_x[k] * osf_517[k];

        t_773[k] = f_8 * nsf_518[k]
                   + f_3 * pc_x[k] * osf_518[k];
    }

#pragma omp simd aligned(t_774, t_775, t_776, pc_x, pc_y, pc_z, nsf_416, nsf_426, nsf_519, \
                         osd0_309, osd1_309, osf_516, osf_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_774[k] = f_8 * nsf_519[k]
                   + f_3 * pc_x[k] * osf_519[k];

        t_775[k] = f_14 * nsf_426[k]
                   + f_1 * osd0_309[k]
                   - f_2 * osd1_309[k]
                   + f_3 * pc_y[k] * osf_516[k];

        t_776[k] = f_17 * nsf_416[k]
                   + f_3 * pc_z[k] * osf_516[k];
    }

#pragma omp simd aligned(t_777, t_778, t_779, pc_y, pc_z, nsf_419, nsf_428, nsf_429, osd0_311, \
                         osd1_311, osf_518, osf_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_777[k] = f_14 * nsf_428[k]
                   + f_4 * osd0_311[k]
                   - f_5 * osd1_311[k]
                   + f_3 * pc_y[k] * osf_518[k];

        t_778[k] = f_14 * nsf_429[k]
                   + f_3 * pc_y[k] * osf_519[k];

        t_779[k] = f_17 * nsf_419[k]
                   + f_1 * osd0_311[k]
                   - f_2 * osd1_311[k]
                   + f_3 * pc_z[k] * osf_519[k];
    }

#pragma omp simd aligned(t_780, t_781, t_782, pc_x, pc_y, pc_z, nsf_420, nsf_430, nsf_520, \
                         osd0_312, osd1_312, osf_520 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_780[k] = f_8 * nsf_520[k]
                   + f_1 * osd0_312[k]
                   - f_2 * osd1_312[k]
                   + f_3 * pc_x[k] * osf_520[k];

        t_781[k] = f_8 * nsf_430[k]
                   + f_3 * pc_y[k] * osf_520[k];

        t_782[k] = f_15 * nsf_420[k]
                   + f_3 * pc_z[k] * osf_520[k];
    }

#pragma omp simd aligned(t_783, t_784, t_785, pc_x, pc_y, nsf_432, nsf_523, nsf_525, osd0_315, \
                         osd0_317, osd1_315, osd1_317, osf_522, osf_523, \
                         osf_525 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_783[k] = f_8 * nsf_523[k]
                   + f_4 * osd0_315[k]
                   - f_5 * osd1_315[k]
                   + f_3 * pc_x[k] * osf_523[k];

        t_784[k] = f_8 * nsf_432[k]
                   + f_3 * pc_y[k] * osf_522[k];

        t_785[k] = f_8 * nsf_525[k]
                   + f_4 * osd0_317[k]
                   - f_5 * osd1_317[k]
                   + f_3 * pc_x[k] * osf_525[k];
    }

#pragma omp simd aligned(t_786, t_787, t_788, t_789, pc_x, nsf_526, nsf_527, nsf_528, nsf_529, \
                         osf_526, osf_527, osf_528, osf_529 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_786[k] = f_8 * nsf_526[k]
                   + f_3 * pc_x[k] * osf_526[k];

        t_787[k] = f_8 * nsf_527[k]
                   + f_3 * pc_x[k] * osf_527[k];

        t_788[k] = f_8 * nsf_528[k]
                   + f_3 * pc_x[k] * osf_528[k];

        t_789[k] = f_8 * nsf_529[k]
                   + f_3 * pc_x[k] * osf_529[k];
    }

#pragma omp simd aligned(t_790, t_791, t_792, pc_y, pc_z, nsf_426, nsf_436, nsf_438, osd0_315, \
                         osd0_317, osd1_315, osd1_317, osf_526, \
                         osf_528 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_790[k] = f_8 * nsf_436[k]
                   + f_1 * osd0_315[k]
                   - f_2 * osd1_315[k]
                   + f_3 * pc_y[k] * osf_526[k];

        t_791[k] = f_15 * nsf_426[k]
                   + f_3 * pc_z[k] * osf_526[k];

        t_792[k] = f_8 * nsf_438[k]
                   + f_4 * osd0_317[k]
                   - f_5 * osd1_317[k]
                   + f_3 * pc_y[k] * osf_528[k];
    }

#pragma omp simd aligned(t_793, t_794, t_795, t_796, pa_y, pc_y, pc_z, nsg0_660, nsf_429, \
                         nsf_439, nsf_440, nsg1_660, osd0_317, osd1_317, osf_529, \
                         osf_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_793[k] = f_8 * nsf_439[k]
                   + f_3 * pc_y[k] * osf_529[k];

        t_794[k] = f_15 * nsf_429[k]
                   + f_1 * osd0_317[k]
                   - f_2 * osd1_317[k]
                   + f_3 * pc_z[k] * osf_529[k];

        t_795[k] = pa_y[k] * nsg0_660[k]
                   - f_6 * pc_y[k] * nsg1_660[k];

        t_796[k] = f_7 * nsf_440[k]
                   + f_3 * pc_y[k] * osf_530[k];
    }

#pragma omp simd aligned(t_797, t_798, t_799, t_800, pa_y, pc_y, pc_z, nsg0_663, nsg0_665, \
                         nsf_430, nsf_441, nsf_442, nsg1_663, nsg1_665, osf_530, \
                         osf_532 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_797[k] = f_13 * nsf_430[k]
                   + f_3 * pc_z[k] * osf_530[k];

        t_798[k] = pa_y[k] * nsg0_663[k]
                   + f_8 * nsf_441[k]
                   - f_6 * pc_y[k] * nsg1_663[k];

        t_799[k] = f_7 * nsf_442[k]
                   + f_3 * pc_y[k] * osf_532[k];

        t_800[k] = pa_y[k] * nsg0_665[k]
                   - f_6 * pc_y[k] * nsg1_665[k];
    }

#pragma omp simd aligned(t_801, t_802, t_803, t_804, pc_x, nsf_536, nsf_537, nsf_538, nsf_539, \
                         osf_536, osf_537, osf_538, osf_539 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_801[k] = f_8 * nsf_536[k]
                   + f_3 * pc_x[k] * osf_536[k];

        t_802[k] = f_8 * nsf_537[k]
                   + f_3 * pc_x[k] * osf_537[k];

        t_803[k] = f_8 * nsf_538[k]
                   + f_3 * pc_x[k] * osf_538[k];

        t_804[k] = f_8 * nsf_539[k]
                   + f_3 * pc_x[k] * osf_539[k];
    }

#pragma omp simd aligned(t_805, t_806, t_807, pc_y, pc_z, nsf_436, nsf_446, nsf_448, osd0_321, \
                         osd0_323, osd1_321, osd1_323, osf_536, \
                         osf_538 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_805[k] = f_7 * nsf_446[k]
                   + f_1 * osd0_321[k]
                   - f_2 * osd1_321[k]
                   + f_3 * pc_y[k] * osf_536[k];

        t_806[k] = f_13 * nsf_436[k]
                   + f_3 * pc_z[k] * osf_536[k];

        t_807[k] = f_7 * nsf_448[k]
                   + f_4 * osd0_323[k]
                   - f_5 * osd1_323[k]
                   + f_3 * pc_y[k] * osf_538[k];
    }

#pragma omp simd aligned(t_808, t_809, t_810, t_811, pa_y, pc_x, pc_y, nsg0_674, nsf_449, \
                         nsf_540, nsg1_674, osd0_324, osd1_324, osf_539, \
                         osf_540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_808[k] = f_7 * nsf_449[k]
                   + f_3 * pc_y[k] * osf_539[k];

        t_809[k] = pa_y[k] * nsg0_674[k]
                   - f_6 * pc_y[k] * nsg1_674[k];

        t_810[k] = f_8 * nsf_540[k]
                   + f_1 * osd0_324[k]
                   - f_2 * osd1_324[k]
                   + f_3 * pc_x[k] * osf_540[k];

        t_811[k] = f_3 * pc_y[k] * osf_540[k];
    }

#pragma omp simd aligned(t_812, t_813, t_814, pc_y, pc_z, nsf_440, osd0_324, osd1_324, \
                         osf_540, osf_541, osf_542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_812[k] = f_12 * nsf_440[k]
                   + f_3 * pc_z[k] * osf_540[k];

        t_813[k] = f_4 * osd0_324[k]
                   - f_5 * osd1_324[k]
                   + f_3 * pc_y[k] * osf_541[k];

        t_814[k] = f_3 * pc_y[k] * osf_542[k];
    }

#pragma omp simd aligned(t_815, t_816, t_817, t_818, pc_x, pc_y, nsf_545, nsf_546, nsf_547, \
                         osd0_329, osd1_329, osf_545, osf_546, \
                         osf_547 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_815[k] = f_8 * nsf_545[k]
                   + f_4 * osd0_329[k]
                   - f_5 * osd1_329[k]
                   + f_3 * pc_x[k] * osf_545[k];

        t_816[k] = f_8 * nsf_546[k]
                   + f_3 * pc_x[k] * osf_546[k];

        t_817[k] = f_8 * nsf_547[k]
                   + f_3 * pc_x[k] * osf_547[k];

        t_818[k] = f_3 * pc_y[k] * osf_545[k];
    }

#pragma omp simd aligned(t_819, t_820, t_821, pc_x, pc_y, nsf_549, osd0_327, osd0_328, \
                         osd1_327, osd1_328, osf_546, osf_547, \
                         osf_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_819[k] = f_8 * nsf_549[k]
                   + f_3 * pc_x[k] * osf_549[k];

        t_820[k] = f_1 * osd0_327[k]
                   - f_2 * osd1_327[k]
                   + f_3 * pc_y[k] * osf_546[k];

        t_821[k] = f_10 * osd0_328[k]
                   - f_11 * osd1_328[k]
                   + f_3 * pc_y[k] * osf_547[k];
    }

#pragma omp simd aligned(t_822, t_823, t_824, t_825, pa_x, pc_x, pc_y, pc_z, nsg0_825, \
                         nsf_449, nsf_550, nsg1_825, osd0_329, osd1_329, osf_548, \
                         osf_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_822[k] = f_4 * osd0_329[k]
                   - f_5 * osd1_329[k]
                   + f_3 * pc_y[k] * osf_548[k];

        t_823[k] = f_3 * pc_y[k] * osf_549[k];

        t_824[k] = f_12 * nsf_449[k]
                   + f_1 * osd0_329[k]
                   - f_2 * osd1_329[k]
                   + f_3 * pc_z[k] * osf_549[k];

        t_825[k] = pa_x[k] * nsg0_825[k]
                   + f_16 * nsf_550[k]
                   - f_6 * pc_x[k] * nsg1_825[k];
    }

#pragma omp simd aligned(t_826, t_827, t_828, t_829, pa_x, pc_x, pc_y, pc_z, nsg0_828, \
                         nsf_450, nsf_553, nsg1_828, osf_550, osf_551 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_826[k] = f_9 * nsf_450[k]
                   + f_3 * pc_y[k] * osf_550[k];

        t_827[k] = f_3 * pc_z[k] * osf_550[k];

        t_828[k] = pa_x[k] * nsg0_828[k]
                   + f_8 * nsf_553[k]
                   - f_6 * pc_x[k] * nsg1_828[k];

        t_829[k] = f_3 * pc_z[k] * osf_551[k];
    }

#pragma omp simd aligned(t_830, t_831, t_832, t_833, pc_x, pc_z, nsf_556, nsf_558, osd0_330, \
                         osd1_330, osf_552, osf_553, osf_556, osf_558 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_830[k] = f_4 * osd0_330[k]
                   - f_5 * osd1_330[k]
                   + f_3 * pc_z[k] * osf_552[k];

        t_831[k] = f_7 * nsf_556[k]
                   + f_3 * pc_x[k] * osf_556[k];

        t_832[k] = f_3 * pc_z[k] * osf_553[k];

        t_833[k] = f_7 * nsf_558[k]
                   + f_3 * pc_x[k] * osf_558[k];
    }

#pragma omp simd aligned(t_834, t_835, t_836, t_837, pa_x, pc_x, pc_z, nsg0_835, nsg0_837, \
                         nsf_559, nsg1_835, nsg1_837, osf_556, \
                         osf_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_834[k] = f_7 * nsf_559[k]
                   + f_3 * pc_x[k] * osf_559[k];

        t_835[k] = pa_x[k] * nsg0_835[k]
                   - f_6 * pc_x[k] * nsg1_835[k];

        t_836[k] = f_3 * pc_z[k] * osf_556[k];

        t_837[k] = pa_x[k] * nsg0_837[k]
                   - f_6 * pc_x[k] * nsg1_837[k];
    }
}

static auto
compute_prim_osg_three_center_electron_repulsion_0_piece7(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsg0,
                                                          const size_t nsf, const size_t nsg1,
                                                          const size_t osf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_3 = p / q;
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 5.0 / q;
    const auto f_12 = 4.5 / q;
    const auto f_13 = 4.0 / q;
    const auto f_14 = 1.5 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.5 / q;

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
    auto *t_945 = buffer.data(target + 945);
    auto *t_946 = buffer.data(target + 946);
    auto *t_947 = buffer.data(target + 947);
    auto *t_948 = buffer.data(target + 948);
    auto *t_949 = buffer.data(target + 949);
    auto *t_950 = buffer.data(target + 950);
    auto *t_951 = buffer.data(target + 951);
    auto *t_952 = buffer.data(target + 952);
    auto *t_953 = buffer.data(target + 953);
    auto *t_954 = buffer.data(target + 954);
    auto *t_955 = buffer.data(target + 955);
    auto *t_956 = buffer.data(target + 956);
    auto *t_957 = buffer.data(target + 957);
    auto *t_958 = buffer.data(target + 958);
    auto *t_959 = buffer.data(target + 959);
    auto *t_960 = buffer.data(target + 960);
    auto *t_961 = buffer.data(target + 961);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsg0_675 = buffer.data(nsg0 + 675);
    const auto *nsg0_678 = buffer.data(nsg0 + 678);
    const auto *nsg0_810 = buffer.data(nsg0 + 810);
    const auto *nsg0_839 = buffer.data(nsg0 + 839);
    const auto *nsg0_845 = buffer.data(nsg0 + 845);
    const auto *nsg0_850 = buffer.data(nsg0 + 850);
    const auto *nsg0_852 = buffer.data(nsg0 + 852);
    const auto *nsg0_854 = buffer.data(nsg0 + 854);
    const auto *nsg0_855 = buffer.data(nsg0 + 855);
    const auto *nsg0_858 = buffer.data(nsg0 + 858);
    const auto *nsg0_860 = buffer.data(nsg0 + 860);
    const auto *nsg0_865 = buffer.data(nsg0 + 865);
    const auto *nsg0_867 = buffer.data(nsg0 + 867);
    const auto *nsg0_869 = buffer.data(nsg0 + 869);
    const auto *nsg0_870 = buffer.data(nsg0 + 870);
    const auto *nsg0_873 = buffer.data(nsg0 + 873);
    const auto *nsg0_875 = buffer.data(nsg0 + 875);
    const auto *nsg0_880 = buffer.data(nsg0 + 880);
    const auto *nsg0_882 = buffer.data(nsg0 + 882);
    const auto *nsg0_884 = buffer.data(nsg0 + 884);
    const auto *nsg0_885 = buffer.data(nsg0 + 885);
    const auto *nsg0_888 = buffer.data(nsg0 + 888);
    const auto *nsg0_890 = buffer.data(nsg0 + 890);
    const auto *nsg0_895 = buffer.data(nsg0 + 895);
    const auto *nsg0_897 = buffer.data(nsg0 + 897);
    const auto *nsg0_899 = buffer.data(nsg0 + 899);
    const auto *nsg0_900 = buffer.data(nsg0 + 900);
    const auto *nsg0_903 = buffer.data(nsg0 + 903);
    const auto *nsg0_905 = buffer.data(nsg0 + 905);
    const auto *nsg0_910 = buffer.data(nsg0 + 910);
    const auto *nsg0_912 = buffer.data(nsg0 + 912);
    const auto *nsg0_914 = buffer.data(nsg0 + 914);
    const auto *nsg0_915 = buffer.data(nsg0 + 915);
    const auto *nsg0_918 = buffer.data(nsg0 + 918);
    const auto *nsg0_920 = buffer.data(nsg0 + 920);
    const auto *nsg0_925 = buffer.data(nsg0 + 925);
    const auto *nsg0_927 = buffer.data(nsg0 + 927);
    const auto *nsg0_929 = buffer.data(nsg0 + 929);
    const auto *nsg0_930 = buffer.data(nsg0 + 930);
    const auto *nsg0_933 = buffer.data(nsg0 + 933);
    const auto *nsg0_935 = buffer.data(nsg0 + 935);
    const auto *nsg0_940 = buffer.data(nsg0 + 940);
    const auto *nsg0_942 = buffer.data(nsg0 + 942);
    const auto *nsg0_944 = buffer.data(nsg0 + 944);
    const auto *nsg0_945 = buffer.data(nsg0 + 945);
    const auto *nsg0_948 = buffer.data(nsg0 + 948);
    const auto *nsg0_950 = buffer.data(nsg0 + 950);
    const auto *nsg0_955 = buffer.data(nsg0 + 955);
    const auto *nsg0_957 = buffer.data(nsg0 + 957);
    const auto *nsg0_959 = buffer.data(nsg0 + 959);

    const auto *nsf_450 = buffer.data(nsf + 450);
    const auto *nsf_456 = buffer.data(nsf + 456);
    const auto *nsf_459 = buffer.data(nsf + 459);
    const auto *nsf_460 = buffer.data(nsf + 460);
    const auto *nsf_462 = buffer.data(nsf + 462);
    const auto *nsf_466 = buffer.data(nsf + 466);
    const auto *nsf_469 = buffer.data(nsf + 469);
    const auto *nsf_470 = buffer.data(nsf + 470);
    const auto *nsf_472 = buffer.data(nsf + 472);
    const auto *nsf_476 = buffer.data(nsf + 476);
    const auto *nsf_479 = buffer.data(nsf + 479);
    const auto *nsf_480 = buffer.data(nsf + 480);
    const auto *nsf_482 = buffer.data(nsf + 482);
    const auto *nsf_486 = buffer.data(nsf + 486);
    const auto *nsf_489 = buffer.data(nsf + 489);
    const auto *nsf_490 = buffer.data(nsf + 490);
    const auto *nsf_492 = buffer.data(nsf + 492);
    const auto *nsf_496 = buffer.data(nsf + 496);
    const auto *nsf_499 = buffer.data(nsf + 499);
    const auto *nsf_500 = buffer.data(nsf + 500);
    const auto *nsf_502 = buffer.data(nsf + 502);
    const auto *nsf_506 = buffer.data(nsf + 506);
    const auto *nsf_509 = buffer.data(nsf + 509);
    const auto *nsf_510 = buffer.data(nsf + 510);
    const auto *nsf_512 = buffer.data(nsf + 512);
    const auto *nsf_516 = buffer.data(nsf + 516);
    const auto *nsf_519 = buffer.data(nsf + 519);
    const auto *nsf_520 = buffer.data(nsf + 520);
    const auto *nsf_522 = buffer.data(nsf + 522);
    const auto *nsf_526 = buffer.data(nsf + 526);
    const auto *nsf_529 = buffer.data(nsf + 529);
    const auto *nsf_530 = buffer.data(nsf + 530);
    const auto *nsf_532 = buffer.data(nsf + 532);
    const auto *nsf_539 = buffer.data(nsf + 539);
    const auto *nsf_540 = buffer.data(nsf + 540);
    const auto *nsf_565 = buffer.data(nsf + 565);
    const auto *nsf_566 = buffer.data(nsf + 566);
    const auto *nsf_567 = buffer.data(nsf + 567);
    const auto *nsf_568 = buffer.data(nsf + 568);
    const auto *nsf_569 = buffer.data(nsf + 569);
    const auto *nsf_570 = buffer.data(nsf + 570);
    const auto *nsf_573 = buffer.data(nsf + 573);
    const auto *nsf_575 = buffer.data(nsf + 575);
    const auto *nsf_576 = buffer.data(nsf + 576);
    const auto *nsf_577 = buffer.data(nsf + 577);
    const auto *nsf_578 = buffer.data(nsf + 578);
    const auto *nsf_579 = buffer.data(nsf + 579);
    const auto *nsf_580 = buffer.data(nsf + 580);
    const auto *nsf_583 = buffer.data(nsf + 583);
    const auto *nsf_585 = buffer.data(nsf + 585);
    const auto *nsf_586 = buffer.data(nsf + 586);
    const auto *nsf_587 = buffer.data(nsf + 587);
    const auto *nsf_588 = buffer.data(nsf + 588);
    const auto *nsf_589 = buffer.data(nsf + 589);
    const auto *nsf_590 = buffer.data(nsf + 590);
    const auto *nsf_593 = buffer.data(nsf + 593);
    const auto *nsf_595 = buffer.data(nsf + 595);
    const auto *nsf_596 = buffer.data(nsf + 596);
    const auto *nsf_597 = buffer.data(nsf + 597);
    const auto *nsf_598 = buffer.data(nsf + 598);
    const auto *nsf_599 = buffer.data(nsf + 599);
    const auto *nsf_600 = buffer.data(nsf + 600);
    const auto *nsf_603 = buffer.data(nsf + 603);
    const auto *nsf_605 = buffer.data(nsf + 605);
    const auto *nsf_606 = buffer.data(nsf + 606);
    const auto *nsf_607 = buffer.data(nsf + 607);
    const auto *nsf_608 = buffer.data(nsf + 608);
    const auto *nsf_609 = buffer.data(nsf + 609);
    const auto *nsf_610 = buffer.data(nsf + 610);
    const auto *nsf_613 = buffer.data(nsf + 613);
    const auto *nsf_615 = buffer.data(nsf + 615);
    const auto *nsf_616 = buffer.data(nsf + 616);
    const auto *nsf_617 = buffer.data(nsf + 617);
    const auto *nsf_618 = buffer.data(nsf + 618);
    const auto *nsf_619 = buffer.data(nsf + 619);
    const auto *nsf_620 = buffer.data(nsf + 620);
    const auto *nsf_623 = buffer.data(nsf + 623);
    const auto *nsf_625 = buffer.data(nsf + 625);
    const auto *nsf_626 = buffer.data(nsf + 626);
    const auto *nsf_627 = buffer.data(nsf + 627);
    const auto *nsf_628 = buffer.data(nsf + 628);
    const auto *nsf_629 = buffer.data(nsf + 629);
    const auto *nsf_630 = buffer.data(nsf + 630);
    const auto *nsf_633 = buffer.data(nsf + 633);
    const auto *nsf_635 = buffer.data(nsf + 635);
    const auto *nsf_636 = buffer.data(nsf + 636);
    const auto *nsf_637 = buffer.data(nsf + 637);
    const auto *nsf_638 = buffer.data(nsf + 638);
    const auto *nsf_639 = buffer.data(nsf + 639);

    const auto *nsg1_675 = buffer.data(nsg1 + 675);
    const auto *nsg1_678 = buffer.data(nsg1 + 678);
    const auto *nsg1_810 = buffer.data(nsg1 + 810);
    const auto *nsg1_839 = buffer.data(nsg1 + 839);
    const auto *nsg1_845 = buffer.data(nsg1 + 845);
    const auto *nsg1_850 = buffer.data(nsg1 + 850);
    const auto *nsg1_852 = buffer.data(nsg1 + 852);
    const auto *nsg1_854 = buffer.data(nsg1 + 854);
    const auto *nsg1_855 = buffer.data(nsg1 + 855);
    const auto *nsg1_858 = buffer.data(nsg1 + 858);
    const auto *nsg1_860 = buffer.data(nsg1 + 860);
    const auto *nsg1_865 = buffer.data(nsg1 + 865);
    const auto *nsg1_867 = buffer.data(nsg1 + 867);
    const auto *nsg1_869 = buffer.data(nsg1 + 869);
    const auto *nsg1_870 = buffer.data(nsg1 + 870);
    const auto *nsg1_873 = buffer.data(nsg1 + 873);
    const auto *nsg1_875 = buffer.data(nsg1 + 875);
    const auto *nsg1_880 = buffer.data(nsg1 + 880);
    const auto *nsg1_882 = buffer.data(nsg1 + 882);
    const auto *nsg1_884 = buffer.data(nsg1 + 884);
    const auto *nsg1_885 = buffer.data(nsg1 + 885);
    const auto *nsg1_888 = buffer.data(nsg1 + 888);
    const auto *nsg1_890 = buffer.data(nsg1 + 890);
    const auto *nsg1_895 = buffer.data(nsg1 + 895);
    const auto *nsg1_897 = buffer.data(nsg1 + 897);
    const auto *nsg1_899 = buffer.data(nsg1 + 899);
    const auto *nsg1_900 = buffer.data(nsg1 + 900);
    const auto *nsg1_903 = buffer.data(nsg1 + 903);
    const auto *nsg1_905 = buffer.data(nsg1 + 905);
    const auto *nsg1_910 = buffer.data(nsg1 + 910);
    const auto *nsg1_912 = buffer.data(nsg1 + 912);
    const auto *nsg1_914 = buffer.data(nsg1 + 914);
    const auto *nsg1_915 = buffer.data(nsg1 + 915);
    const auto *nsg1_918 = buffer.data(nsg1 + 918);
    const auto *nsg1_920 = buffer.data(nsg1 + 920);
    const auto *nsg1_925 = buffer.data(nsg1 + 925);
    const auto *nsg1_927 = buffer.data(nsg1 + 927);
    const auto *nsg1_929 = buffer.data(nsg1 + 929);
    const auto *nsg1_930 = buffer.data(nsg1 + 930);
    const auto *nsg1_933 = buffer.data(nsg1 + 933);
    const auto *nsg1_935 = buffer.data(nsg1 + 935);
    const auto *nsg1_940 = buffer.data(nsg1 + 940);
    const auto *nsg1_942 = buffer.data(nsg1 + 942);
    const auto *nsg1_944 = buffer.data(nsg1 + 944);
    const auto *nsg1_945 = buffer.data(nsg1 + 945);
    const auto *nsg1_948 = buffer.data(nsg1 + 948);
    const auto *nsg1_950 = buffer.data(nsg1 + 950);
    const auto *nsg1_955 = buffer.data(nsg1 + 955);
    const auto *nsg1_957 = buffer.data(nsg1 + 957);
    const auto *nsg1_959 = buffer.data(nsg1 + 959);

    const auto *osf_559 = buffer.data(osf + 559);
    const auto *osf_560 = buffer.data(osf + 560);
    const auto *osf_562 = buffer.data(osf + 562);
    const auto *osf_566 = buffer.data(osf + 566);
    const auto *osf_567 = buffer.data(osf + 567);
    const auto *osf_568 = buffer.data(osf + 568);
    const auto *osf_569 = buffer.data(osf + 569);
    const auto *osf_570 = buffer.data(osf + 570);
    const auto *osf_572 = buffer.data(osf + 572);
    const auto *osf_576 = buffer.data(osf + 576);
    const auto *osf_577 = buffer.data(osf + 577);
    const auto *osf_578 = buffer.data(osf + 578);
    const auto *osf_579 = buffer.data(osf + 579);
    const auto *osf_580 = buffer.data(osf + 580);
    const auto *osf_582 = buffer.data(osf + 582);
    const auto *osf_586 = buffer.data(osf + 586);
    const auto *osf_587 = buffer.data(osf + 587);
    const auto *osf_588 = buffer.data(osf + 588);
    const auto *osf_589 = buffer.data(osf + 589);
    const auto *osf_590 = buffer.data(osf + 590);
    const auto *osf_592 = buffer.data(osf + 592);
    const auto *osf_596 = buffer.data(osf + 596);
    const auto *osf_597 = buffer.data(osf + 597);
    const auto *osf_598 = buffer.data(osf + 598);
    const auto *osf_599 = buffer.data(osf + 599);
    const auto *osf_600 = buffer.data(osf + 600);
    const auto *osf_602 = buffer.data(osf + 602);
    const auto *osf_606 = buffer.data(osf + 606);
    const auto *osf_607 = buffer.data(osf + 607);
    const auto *osf_608 = buffer.data(osf + 608);
    const auto *osf_609 = buffer.data(osf + 609);
    const auto *osf_610 = buffer.data(osf + 610);
    const auto *osf_612 = buffer.data(osf + 612);
    const auto *osf_616 = buffer.data(osf + 616);
    const auto *osf_617 = buffer.data(osf + 617);
    const auto *osf_618 = buffer.data(osf + 618);
    const auto *osf_619 = buffer.data(osf + 619);
    const auto *osf_620 = buffer.data(osf + 620);
    const auto *osf_622 = buffer.data(osf + 622);
    const auto *osf_626 = buffer.data(osf + 626);
    const auto *osf_627 = buffer.data(osf + 627);
    const auto *osf_628 = buffer.data(osf + 628);
    const auto *osf_629 = buffer.data(osf + 629);
    const auto *osf_630 = buffer.data(osf + 630);
    const auto *osf_632 = buffer.data(osf + 632);
    const auto *osf_636 = buffer.data(osf + 636);
    const auto *osf_637 = buffer.data(osf + 637);
    const auto *osf_638 = buffer.data(osf + 638);
    const auto *osf_639 = buffer.data(osf + 639);
    const auto *osf_640 = buffer.data(osf + 640);

#pragma omp simd aligned(t_838, t_839, t_840, pa_x, pa_z, pc_x, pc_y, pc_z, nsg0_675, \
                         nsg0_839, nsf_459, nsg1_675, nsg1_839, \
                         osf_559 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_838[k] = f_9 * nsf_459[k]
                   + f_3 * pc_y[k] * osf_559[k];

        t_839[k] = pa_x[k] * nsg0_839[k]
                   - f_6 * pc_x[k] * nsg1_839[k];

        t_840[k] = pa_z[k] * nsg0_675[k]
                   - f_6 * pc_z[k] * nsg1_675[k];
    }

#pragma omp simd aligned(t_841, t_842, t_843, t_844, pa_z, pc_y, pc_z, nsg0_678, nsf_450, \
                         nsf_460, nsf_462, nsg1_678, osf_560, osf_562 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_841[k] = f_12 * nsf_460[k]
                   + f_3 * pc_y[k] * osf_560[k];

        t_842[k] = f_7 * nsf_450[k]
                   + f_3 * pc_z[k] * osf_560[k];

        t_843[k] = pa_z[k] * nsg0_678[k]
                   - f_6 * pc_z[k] * nsg1_678[k];

        t_844[k] = f_12 * nsf_462[k]
                   + f_3 * pc_y[k] * osf_562[k];
    }

#pragma omp simd aligned(t_845, t_846, t_847, t_848, pa_x, pc_x, nsg0_845, nsf_565, nsf_566, \
                         nsf_567, nsf_568, nsg1_845, osf_566, osf_567, \
                         osf_568 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_845[k] = pa_x[k] * nsg0_845[k]
                   + f_8 * nsf_565[k]
                   - f_6 * pc_x[k] * nsg1_845[k];

        t_846[k] = f_7 * nsf_566[k]
                   + f_3 * pc_x[k] * osf_566[k];

        t_847[k] = f_7 * nsf_567[k]
                   + f_3 * pc_x[k] * osf_567[k];

        t_848[k] = f_7 * nsf_568[k]
                   + f_3 * pc_x[k] * osf_568[k];
    }

#pragma omp simd aligned(t_849, t_850, t_851, t_852, pa_x, pc_x, pc_z, nsg0_850, nsg0_852, \
                         nsf_456, nsf_569, nsg1_850, nsg1_852, osf_566, \
                         osf_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_849[k] = f_7 * nsf_569[k]
                   + f_3 * pc_x[k] * osf_569[k];

        t_850[k] = pa_x[k] * nsg0_850[k]
                   - f_6 * pc_x[k] * nsg1_850[k];

        t_851[k] = f_7 * nsf_456[k]
                   + f_3 * pc_z[k] * osf_566[k];

        t_852[k] = pa_x[k] * nsg0_852[k]
                   - f_6 * pc_x[k] * nsg1_852[k];
    }

#pragma omp simd aligned(t_853, t_854, t_855, t_856, pa_x, pc_x, pc_y, nsg0_854, nsg0_855, \
                         nsf_469, nsf_470, nsf_570, nsg1_854, nsg1_855, osf_569, \
                         osf_570 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_853[k] = f_12 * nsf_469[k]
                   + f_3 * pc_y[k] * osf_569[k];

        t_854[k] = pa_x[k] * nsg0_854[k]
                   - f_6 * pc_x[k] * nsg1_854[k];

        t_855[k] = pa_x[k] * nsg0_855[k]
                   + f_16 * nsf_570[k]
                   - f_6 * pc_x[k] * nsg1_855[k];

        t_856[k] = f_13 * nsf_470[k]
                   + f_3 * pc_y[k] * osf_570[k];
    }

#pragma omp simd aligned(t_857, t_858, t_859, pa_x, pc_x, pc_y, pc_z, nsg0_858, nsf_460, \
                         nsf_472, nsf_573, nsg1_858, osf_570, osf_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_857[k] = f_8 * nsf_460[k]
                   + f_3 * pc_z[k] * osf_570[k];

        t_858[k] = pa_x[k] * nsg0_858[k]
                   + f_8 * nsf_573[k]
                   - f_6 * pc_x[k] * nsg1_858[k];

        t_859[k] = f_13 * nsf_472[k]
                   + f_3 * pc_y[k] * osf_572[k];
    }

#pragma omp simd aligned(t_860, t_861, t_862, t_863, pa_x, pc_x, nsg0_860, nsf_575, nsf_576, \
                         nsf_577, nsf_578, nsg1_860, osf_576, osf_577, \
                         osf_578 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_860[k] = pa_x[k] * nsg0_860[k]
                   + f_8 * nsf_575[k]
                   - f_6 * pc_x[k] * nsg1_860[k];

        t_861[k] = f_7 * nsf_576[k]
                   + f_3 * pc_x[k] * osf_576[k];

        t_862[k] = f_7 * nsf_577[k]
                   + f_3 * pc_x[k] * osf_577[k];

        t_863[k] = f_7 * nsf_578[k]
                   + f_3 * pc_x[k] * osf_578[k];
    }

#pragma omp simd aligned(t_864, t_865, t_866, t_867, pa_x, pc_x, pc_z, nsg0_865, nsg0_867, \
                         nsf_466, nsf_579, nsg1_865, nsg1_867, osf_576, \
                         osf_579 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_864[k] = f_7 * nsf_579[k]
                   + f_3 * pc_x[k] * osf_579[k];

        t_865[k] = pa_x[k] * nsg0_865[k]
                   - f_6 * pc_x[k] * nsg1_865[k];

        t_866[k] = f_8 * nsf_466[k]
                   + f_3 * pc_z[k] * osf_576[k];

        t_867[k] = pa_x[k] * nsg0_867[k]
                   - f_6 * pc_x[k] * nsg1_867[k];
    }

#pragma omp simd aligned(t_868, t_869, t_870, t_871, pa_x, pc_x, pc_y, nsg0_869, nsg0_870, \
                         nsf_479, nsf_480, nsf_580, nsg1_869, nsg1_870, osf_579, \
                         osf_580 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_868[k] = f_13 * nsf_479[k]
                   + f_3 * pc_y[k] * osf_579[k];

        t_869[k] = pa_x[k] * nsg0_869[k]
                   - f_6 * pc_x[k] * nsg1_869[k];

        t_870[k] = pa_x[k] * nsg0_870[k]
                   + f_16 * nsf_580[k]
                   - f_6 * pc_x[k] * nsg1_870[k];

        t_871[k] = f_15 * nsf_480[k]
                   + f_3 * pc_y[k] * osf_580[k];
    }

#pragma omp simd aligned(t_872, t_873, t_874, pa_x, pc_x, pc_y, pc_z, nsg0_873, nsf_470, \
                         nsf_482, nsf_583, nsg1_873, osf_580, osf_582 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_872[k] = f_14 * nsf_470[k]
                   + f_3 * pc_z[k] * osf_580[k];

        t_873[k] = pa_x[k] * nsg0_873[k]
                   + f_8 * nsf_583[k]
                   - f_6 * pc_x[k] * nsg1_873[k];

        t_874[k] = f_15 * nsf_482[k]
                   + f_3 * pc_y[k] * osf_582[k];
    }

#pragma omp simd aligned(t_875, t_876, t_877, t_878, pa_x, pc_x, nsg0_875, nsf_585, nsf_586, \
                         nsf_587, nsf_588, nsg1_875, osf_586, osf_587, \
                         osf_588 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_875[k] = pa_x[k] * nsg0_875[k]
                   + f_8 * nsf_585[k]
                   - f_6 * pc_x[k] * nsg1_875[k];

        t_876[k] = f_7 * nsf_586[k]
                   + f_3 * pc_x[k] * osf_586[k];

        t_877[k] = f_7 * nsf_587[k]
                   + f_3 * pc_x[k] * osf_587[k];

        t_878[k] = f_7 * nsf_588[k]
                   + f_3 * pc_x[k] * osf_588[k];
    }

#pragma omp simd aligned(t_879, t_880, t_881, t_882, pa_x, pc_x, pc_z, nsg0_880, nsg0_882, \
                         nsf_476, nsf_589, nsg1_880, nsg1_882, osf_586, \
                         osf_589 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_879[k] = f_7 * nsf_589[k]
                   + f_3 * pc_x[k] * osf_589[k];

        t_880[k] = pa_x[k] * nsg0_880[k]
                   - f_6 * pc_x[k] * nsg1_880[k];

        t_881[k] = f_14 * nsf_476[k]
                   + f_3 * pc_z[k] * osf_586[k];

        t_882[k] = pa_x[k] * nsg0_882[k]
                   - f_6 * pc_x[k] * nsg1_882[k];
    }

#pragma omp simd aligned(t_883, t_884, t_885, t_886, pa_x, pc_x, pc_y, nsg0_884, nsg0_885, \
                         nsf_489, nsf_490, nsf_590, nsg1_884, nsg1_885, osf_589, \
                         osf_590 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_883[k] = f_15 * nsf_489[k]
                   + f_3 * pc_y[k] * osf_589[k];

        t_884[k] = pa_x[k] * nsg0_884[k]
                   - f_6 * pc_x[k] * nsg1_884[k];

        t_885[k] = pa_x[k] * nsg0_885[k]
                   + f_16 * nsf_590[k]
                   - f_6 * pc_x[k] * nsg1_885[k];

        t_886[k] = f_17 * nsf_490[k]
                   + f_3 * pc_y[k] * osf_590[k];
    }

#pragma omp simd aligned(t_887, t_888, t_889, pa_x, pc_x, pc_y, pc_z, nsg0_888, nsf_480, \
                         nsf_492, nsf_593, nsg1_888, osf_590, osf_592 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_887[k] = f_16 * nsf_480[k]
                   + f_3 * pc_z[k] * osf_590[k];

        t_888[k] = pa_x[k] * nsg0_888[k]
                   + f_8 * nsf_593[k]
                   - f_6 * pc_x[k] * nsg1_888[k];

        t_889[k] = f_17 * nsf_492[k]
                   + f_3 * pc_y[k] * osf_592[k];
    }

#pragma omp simd aligned(t_890, t_891, t_892, t_893, pa_x, pc_x, nsg0_890, nsf_595, nsf_596, \
                         nsf_597, nsf_598, nsg1_890, osf_596, osf_597, \
                         osf_598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_890[k] = pa_x[k] * nsg0_890[k]
                   + f_8 * nsf_595[k]
                   - f_6 * pc_x[k] * nsg1_890[k];

        t_891[k] = f_7 * nsf_596[k]
                   + f_3 * pc_x[k] * osf_596[k];

        t_892[k] = f_7 * nsf_597[k]
                   + f_3 * pc_x[k] * osf_597[k];

        t_893[k] = f_7 * nsf_598[k]
                   + f_3 * pc_x[k] * osf_598[k];
    }

#pragma omp simd aligned(t_894, t_895, t_896, t_897, pa_x, pc_x, pc_z, nsg0_895, nsg0_897, \
                         nsf_486, nsf_599, nsg1_895, nsg1_897, osf_596, \
                         osf_599 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_894[k] = f_7 * nsf_599[k]
                   + f_3 * pc_x[k] * osf_599[k];

        t_895[k] = pa_x[k] * nsg0_895[k]
                   - f_6 * pc_x[k] * nsg1_895[k];

        t_896[k] = f_16 * nsf_486[k]
                   + f_3 * pc_z[k] * osf_596[k];

        t_897[k] = pa_x[k] * nsg0_897[k]
                   - f_6 * pc_x[k] * nsg1_897[k];
    }

#pragma omp simd aligned(t_898, t_899, t_900, t_901, pa_x, pc_x, pc_y, nsg0_899, nsg0_900, \
                         nsf_499, nsf_500, nsf_600, nsg1_899, nsg1_900, osf_599, \
                         osf_600 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_898[k] = f_17 * nsf_499[k]
                   + f_3 * pc_y[k] * osf_599[k];

        t_899[k] = pa_x[k] * nsg0_899[k]
                   - f_6 * pc_x[k] * nsg1_899[k];

        t_900[k] = pa_x[k] * nsg0_900[k]
                   + f_16 * nsf_600[k]
                   - f_6 * pc_x[k] * nsg1_900[k];

        t_901[k] = f_18 * nsf_500[k]
                   + f_3 * pc_y[k] * osf_600[k];
    }

#pragma omp simd aligned(t_902, t_903, t_904, pa_x, pc_x, pc_y, pc_z, nsg0_903, nsf_490, \
                         nsf_502, nsf_603, nsg1_903, osf_600, osf_602 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_902[k] = f_18 * nsf_490[k]
                   + f_3 * pc_z[k] * osf_600[k];

        t_903[k] = pa_x[k] * nsg0_903[k]
                   + f_8 * nsf_603[k]
                   - f_6 * pc_x[k] * nsg1_903[k];

        t_904[k] = f_18 * nsf_502[k]
                   + f_3 * pc_y[k] * osf_602[k];
    }

#pragma omp simd aligned(t_905, t_906, t_907, t_908, pa_x, pc_x, nsg0_905, nsf_605, nsf_606, \
                         nsf_607, nsf_608, nsg1_905, osf_606, osf_607, \
                         osf_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_905[k] = pa_x[k] * nsg0_905[k]
                   + f_8 * nsf_605[k]
                   - f_6 * pc_x[k] * nsg1_905[k];

        t_906[k] = f_7 * nsf_606[k]
                   + f_3 * pc_x[k] * osf_606[k];

        t_907[k] = f_7 * nsf_607[k]
                   + f_3 * pc_x[k] * osf_607[k];

        t_908[k] = f_7 * nsf_608[k]
                   + f_3 * pc_x[k] * osf_608[k];
    }

#pragma omp simd aligned(t_909, t_910, t_911, t_912, pa_x, pc_x, pc_z, nsg0_910, nsg0_912, \
                         nsf_496, nsf_609, nsg1_910, nsg1_912, osf_606, \
                         osf_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_909[k] = f_7 * nsf_609[k]
                   + f_3 * pc_x[k] * osf_609[k];

        t_910[k] = pa_x[k] * nsg0_910[k]
                   - f_6 * pc_x[k] * nsg1_910[k];

        t_911[k] = f_18 * nsf_496[k]
                   + f_3 * pc_z[k] * osf_606[k];

        t_912[k] = pa_x[k] * nsg0_912[k]
                   - f_6 * pc_x[k] * nsg1_912[k];
    }

#pragma omp simd aligned(t_913, t_914, t_915, t_916, pa_x, pc_x, pc_y, nsg0_914, nsg0_915, \
                         nsf_509, nsf_510, nsf_610, nsg1_914, nsg1_915, osf_609, \
                         osf_610 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_913[k] = f_18 * nsf_509[k]
                   + f_3 * pc_y[k] * osf_609[k];

        t_914[k] = pa_x[k] * nsg0_914[k]
                   - f_6 * pc_x[k] * nsg1_914[k];

        t_915[k] = pa_x[k] * nsg0_915[k]
                   + f_16 * nsf_610[k]
                   - f_6 * pc_x[k] * nsg1_915[k];

        t_916[k] = f_16 * nsf_510[k]
                   + f_3 * pc_y[k] * osf_610[k];
    }

#pragma omp simd aligned(t_917, t_918, t_919, pa_x, pc_x, pc_y, pc_z, nsg0_918, nsf_500, \
                         nsf_512, nsf_613, nsg1_918, osf_610, osf_612 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_917[k] = f_17 * nsf_500[k]
                   + f_3 * pc_z[k] * osf_610[k];

        t_918[k] = pa_x[k] * nsg0_918[k]
                   + f_8 * nsf_613[k]
                   - f_6 * pc_x[k] * nsg1_918[k];

        t_919[k] = f_16 * nsf_512[k]
                   + f_3 * pc_y[k] * osf_612[k];
    }

#pragma omp simd aligned(t_920, t_921, t_922, t_923, pa_x, pc_x, nsg0_920, nsf_615, nsf_616, \
                         nsf_617, nsf_618, nsg1_920, osf_616, osf_617, \
                         osf_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_920[k] = pa_x[k] * nsg0_920[k]
                   + f_8 * nsf_615[k]
                   - f_6 * pc_x[k] * nsg1_920[k];

        t_921[k] = f_7 * nsf_616[k]
                   + f_3 * pc_x[k] * osf_616[k];

        t_922[k] = f_7 * nsf_617[k]
                   + f_3 * pc_x[k] * osf_617[k];

        t_923[k] = f_7 * nsf_618[k]
                   + f_3 * pc_x[k] * osf_618[k];
    }

#pragma omp simd aligned(t_924, t_925, t_926, t_927, pa_x, pc_x, pc_z, nsg0_925, nsg0_927, \
                         nsf_506, nsf_619, nsg1_925, nsg1_927, osf_616, \
                         osf_619 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_924[k] = f_7 * nsf_619[k]
                   + f_3 * pc_x[k] * osf_619[k];

        t_925[k] = pa_x[k] * nsg0_925[k]
                   - f_6 * pc_x[k] * nsg1_925[k];

        t_926[k] = f_17 * nsf_506[k]
                   + f_3 * pc_z[k] * osf_616[k];

        t_927[k] = pa_x[k] * nsg0_927[k]
                   - f_6 * pc_x[k] * nsg1_927[k];
    }

#pragma omp simd aligned(t_928, t_929, t_930, t_931, pa_x, pc_x, pc_y, nsg0_929, nsg0_930, \
                         nsf_519, nsf_520, nsf_620, nsg1_929, nsg1_930, osf_619, \
                         osf_620 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_928[k] = f_16 * nsf_519[k]
                   + f_3 * pc_y[k] * osf_619[k];

        t_929[k] = pa_x[k] * nsg0_929[k]
                   - f_6 * pc_x[k] * nsg1_929[k];

        t_930[k] = pa_x[k] * nsg0_930[k]
                   + f_16 * nsf_620[k]
                   - f_6 * pc_x[k] * nsg1_930[k];

        t_931[k] = f_14 * nsf_520[k]
                   + f_3 * pc_y[k] * osf_620[k];
    }

#pragma omp simd aligned(t_932, t_933, t_934, pa_x, pc_x, pc_y, pc_z, nsg0_933, nsf_510, \
                         nsf_522, nsf_623, nsg1_933, osf_620, osf_622 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_932[k] = f_15 * nsf_510[k]
                   + f_3 * pc_z[k] * osf_620[k];

        t_933[k] = pa_x[k] * nsg0_933[k]
                   + f_8 * nsf_623[k]
                   - f_6 * pc_x[k] * nsg1_933[k];

        t_934[k] = f_14 * nsf_522[k]
                   + f_3 * pc_y[k] * osf_622[k];
    }

#pragma omp simd aligned(t_935, t_936, t_937, t_938, pa_x, pc_x, nsg0_935, nsf_625, nsf_626, \
                         nsf_627, nsf_628, nsg1_935, osf_626, osf_627, \
                         osf_628 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_935[k] = pa_x[k] * nsg0_935[k]
                   + f_8 * nsf_625[k]
                   - f_6 * pc_x[k] * nsg1_935[k];

        t_936[k] = f_7 * nsf_626[k]
                   + f_3 * pc_x[k] * osf_626[k];

        t_937[k] = f_7 * nsf_627[k]
                   + f_3 * pc_x[k] * osf_627[k];

        t_938[k] = f_7 * nsf_628[k]
                   + f_3 * pc_x[k] * osf_628[k];
    }

#pragma omp simd aligned(t_939, t_940, t_941, t_942, pa_x, pc_x, pc_z, nsg0_940, nsg0_942, \
                         nsf_516, nsf_629, nsg1_940, nsg1_942, osf_626, \
                         osf_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_939[k] = f_7 * nsf_629[k]
                   + f_3 * pc_x[k] * osf_629[k];

        t_940[k] = pa_x[k] * nsg0_940[k]
                   - f_6 * pc_x[k] * nsg1_940[k];

        t_941[k] = f_15 * nsf_516[k]
                   + f_3 * pc_z[k] * osf_626[k];

        t_942[k] = pa_x[k] * nsg0_942[k]
                   - f_6 * pc_x[k] * nsg1_942[k];
    }

#pragma omp simd aligned(t_943, t_944, t_945, t_946, pa_x, pc_x, pc_y, nsg0_944, nsg0_945, \
                         nsf_529, nsf_530, nsf_630, nsg1_944, nsg1_945, osf_629, \
                         osf_630 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_943[k] = f_14 * nsf_529[k]
                   + f_3 * pc_y[k] * osf_629[k];

        t_944[k] = pa_x[k] * nsg0_944[k]
                   - f_6 * pc_x[k] * nsg1_944[k];

        t_945[k] = pa_x[k] * nsg0_945[k]
                   + f_16 * nsf_630[k]
                   - f_6 * pc_x[k] * nsg1_945[k];

        t_946[k] = f_8 * nsf_530[k]
                   + f_3 * pc_y[k] * osf_630[k];
    }

#pragma omp simd aligned(t_947, t_948, t_949, pa_x, pc_x, pc_y, pc_z, nsg0_948, nsf_520, \
                         nsf_532, nsf_633, nsg1_948, osf_630, osf_632 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_947[k] = f_13 * nsf_520[k]
                   + f_3 * pc_z[k] * osf_630[k];

        t_948[k] = pa_x[k] * nsg0_948[k]
                   + f_8 * nsf_633[k]
                   - f_6 * pc_x[k] * nsg1_948[k];

        t_949[k] = f_8 * nsf_532[k]
                   + f_3 * pc_y[k] * osf_632[k];
    }

#pragma omp simd aligned(t_950, t_951, t_952, t_953, pa_x, pc_x, nsg0_950, nsf_635, nsf_636, \
                         nsf_637, nsf_638, nsg1_950, osf_636, osf_637, \
                         osf_638 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_950[k] = pa_x[k] * nsg0_950[k]
                   + f_8 * nsf_635[k]
                   - f_6 * pc_x[k] * nsg1_950[k];

        t_951[k] = f_7 * nsf_636[k]
                   + f_3 * pc_x[k] * osf_636[k];

        t_952[k] = f_7 * nsf_637[k]
                   + f_3 * pc_x[k] * osf_637[k];

        t_953[k] = f_7 * nsf_638[k]
                   + f_3 * pc_x[k] * osf_638[k];
    }

#pragma omp simd aligned(t_954, t_955, t_956, t_957, pa_x, pc_x, pc_z, nsg0_955, nsg0_957, \
                         nsf_526, nsf_639, nsg1_955, nsg1_957, osf_636, \
                         osf_639 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_954[k] = f_7 * nsf_639[k]
                   + f_3 * pc_x[k] * osf_639[k];

        t_955[k] = pa_x[k] * nsg0_955[k]
                   - f_6 * pc_x[k] * nsg1_955[k];

        t_956[k] = f_13 * nsf_526[k]
                   + f_3 * pc_z[k] * osf_636[k];

        t_957[k] = pa_x[k] * nsg0_957[k]
                   - f_6 * pc_x[k] * nsg1_957[k];
    }

#pragma omp simd aligned(t_958, t_959, t_960, t_961, pa_x, pa_y, pc_x, pc_y, nsg0_810, \
                         nsg0_959, nsf_539, nsf_540, nsg1_810, nsg1_959, osf_639, \
                         osf_640 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_958[k] = f_8 * nsf_539[k]
                   + f_3 * pc_y[k] * osf_639[k];

        t_959[k] = pa_x[k] * nsg0_959[k]
                   - f_6 * pc_x[k] * nsg1_959[k];

        t_960[k] = pa_y[k] * nsg0_810[k]
                   - f_6 * pc_y[k] * nsg1_810[k];

        t_961[k] = f_7 * nsf_540[k]
                   + f_3 * pc_y[k] * osf_640[k];
    }
}

static auto
compute_prim_osg_three_center_electron_repulsion_0_piece8(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsg0,
                                                          const size_t nsf, const size_t nsg1,
                                                          const size_t osd0, const size_t osd1,
                                                          const size_t osf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 5.0 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 4.5 / q;
    const auto f_13 = 4.0 / q;
    const auto f_14 = 1.5 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.5 / q;

    auto *t_962 = buffer.data(target + 962);
    auto *t_963 = buffer.data(target + 963);
    auto *t_964 = buffer.data(target + 964);
    auto *t_965 = buffer.data(target + 965);
    auto *t_966 = buffer.data(target + 966);
    auto *t_967 = buffer.data(target + 967);
    auto *t_968 = buffer.data(target + 968);
    auto *t_969 = buffer.data(target + 969);
    auto *t_970 = buffer.data(target + 970);
    auto *t_971 = buffer.data(target + 971);
    auto *t_972 = buffer.data(target + 972);
    auto *t_973 = buffer.data(target + 973);
    auto *t_974 = buffer.data(target + 974);
    auto *t_975 = buffer.data(target + 975);
    auto *t_976 = buffer.data(target + 976);
    auto *t_977 = buffer.data(target + 977);
    auto *t_978 = buffer.data(target + 978);
    auto *t_979 = buffer.data(target + 979);
    auto *t_980 = buffer.data(target + 980);
    auto *t_981 = buffer.data(target + 981);
    auto *t_982 = buffer.data(target + 982);
    auto *t_983 = buffer.data(target + 983);
    auto *t_984 = buffer.data(target + 984);
    auto *t_985 = buffer.data(target + 985);
    auto *t_986 = buffer.data(target + 986);
    auto *t_987 = buffer.data(target + 987);
    auto *t_988 = buffer.data(target + 988);
    auto *t_989 = buffer.data(target + 989);
    auto *t_990 = buffer.data(target + 990);
    auto *t_991 = buffer.data(target + 991);
    auto *t_992 = buffer.data(target + 992);
    auto *t_993 = buffer.data(target + 993);
    auto *t_994 = buffer.data(target + 994);
    auto *t_995 = buffer.data(target + 995);
    auto *t_996 = buffer.data(target + 996);
    auto *t_997 = buffer.data(target + 997);
    auto *t_998 = buffer.data(target + 998);
    auto *t_999 = buffer.data(target + 999);
    auto *t_1000 = buffer.data(target + 1000);
    auto *t_1001 = buffer.data(target + 1001);
    auto *t_1002 = buffer.data(target + 1002);
    auto *t_1003 = buffer.data(target + 1003);
    auto *t_1004 = buffer.data(target + 1004);
    auto *t_1005 = buffer.data(target + 1005);
    auto *t_1006 = buffer.data(target + 1006);
    auto *t_1007 = buffer.data(target + 1007);
    auto *t_1008 = buffer.data(target + 1008);
    auto *t_1009 = buffer.data(target + 1009);
    auto *t_1010 = buffer.data(target + 1010);
    auto *t_1011 = buffer.data(target + 1011);
    auto *t_1012 = buffer.data(target + 1012);
    auto *t_1013 = buffer.data(target + 1013);
    auto *t_1014 = buffer.data(target + 1014);
    auto *t_1015 = buffer.data(target + 1015);
    auto *t_1016 = buffer.data(target + 1016);
    auto *t_1017 = buffer.data(target + 1017);
    auto *t_1018 = buffer.data(target + 1018);
    auto *t_1019 = buffer.data(target + 1019);
    auto *t_1020 = buffer.data(target + 1020);
    auto *t_1021 = buffer.data(target + 1021);
    auto *t_1022 = buffer.data(target + 1022);
    auto *t_1023 = buffer.data(target + 1023);
    auto *t_1024 = buffer.data(target + 1024);
    auto *t_1025 = buffer.data(target + 1025);
    auto *t_1026 = buffer.data(target + 1026);
    auto *t_1027 = buffer.data(target + 1027);
    auto *t_1028 = buffer.data(target + 1028);
    auto *t_1029 = buffer.data(target + 1029);
    auto *t_1030 = buffer.data(target + 1030);
    auto *t_1031 = buffer.data(target + 1031);
    auto *t_1032 = buffer.data(target + 1032);
    auto *t_1033 = buffer.data(target + 1033);
    auto *t_1034 = buffer.data(target + 1034);
    auto *t_1035 = buffer.data(target + 1035);
    auto *t_1036 = buffer.data(target + 1036);
    auto *t_1037 = buffer.data(target + 1037);
    auto *t_1038 = buffer.data(target + 1038);
    auto *t_1039 = buffer.data(target + 1039);
    auto *t_1040 = buffer.data(target + 1040);
    auto *t_1041 = buffer.data(target + 1041);
    auto *t_1042 = buffer.data(target + 1042);
    auto *t_1043 = buffer.data(target + 1043);
    auto *t_1044 = buffer.data(target + 1044);
    auto *t_1045 = buffer.data(target + 1045);
    auto *t_1046 = buffer.data(target + 1046);
    auto *t_1047 = buffer.data(target + 1047);
    auto *t_1048 = buffer.data(target + 1048);
    auto *t_1049 = buffer.data(target + 1049);
    auto *t_1050 = buffer.data(target + 1050);
    auto *t_1051 = buffer.data(target + 1051);
    auto *t_1052 = buffer.data(target + 1052);
    auto *t_1053 = buffer.data(target + 1053);
    auto *t_1054 = buffer.data(target + 1054);
    auto *t_1055 = buffer.data(target + 1055);
    auto *t_1056 = buffer.data(target + 1056);
    auto *t_1057 = buffer.data(target + 1057);
    auto *t_1058 = buffer.data(target + 1058);
    auto *t_1059 = buffer.data(target + 1059);
    auto *t_1060 = buffer.data(target + 1060);
    auto *t_1061 = buffer.data(target + 1061);
    auto *t_1062 = buffer.data(target + 1062);
    auto *t_1063 = buffer.data(target + 1063);
    auto *t_1064 = buffer.data(target + 1064);
    auto *t_1065 = buffer.data(target + 1065);
    auto *t_1066 = buffer.data(target + 1066);
    auto *t_1067 = buffer.data(target + 1067);
    auto *t_1068 = buffer.data(target + 1068);
    auto *t_1069 = buffer.data(target + 1069);
    auto *t_1070 = buffer.data(target + 1070);
    auto *t_1071 = buffer.data(target + 1071);
    auto *t_1072 = buffer.data(target + 1072);
    auto *t_1073 = buffer.data(target + 1073);
    auto *t_1074 = buffer.data(target + 1074);
    auto *t_1075 = buffer.data(target + 1075);
    auto *t_1076 = buffer.data(target + 1076);
    auto *t_1077 = buffer.data(target + 1077);
    auto *t_1078 = buffer.data(target + 1078);
    auto *t_1079 = buffer.data(target + 1079);
    auto *t_1080 = buffer.data(target + 1080);
    auto *t_1081 = buffer.data(target + 1081);
    auto *t_1082 = buffer.data(target + 1082);
    auto *t_1083 = buffer.data(target + 1083);
    auto *t_1084 = buffer.data(target + 1084);
    auto *t_1085 = buffer.data(target + 1085);
    auto *t_1086 = buffer.data(target + 1086);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsg0_815 = buffer.data(nsg0 + 815);
    const auto *nsg0_825 = buffer.data(nsg0 + 825);
    const auto *nsg0_826 = buffer.data(nsg0 + 826);
    const auto *nsg0_828 = buffer.data(nsg0 + 828);
    const auto *nsg0_835 = buffer.data(nsg0 + 835);
    const auto *nsg0_837 = buffer.data(nsg0 + 837);
    const auto *nsg0_963 = buffer.data(nsg0 + 963);
    const auto *nsg0_970 = buffer.data(nsg0 + 970);
    const auto *nsg0_972 = buffer.data(nsg0 + 972);
    const auto *nsg0_974 = buffer.data(nsg0 + 974);
    const auto *nsg0_975 = buffer.data(nsg0 + 975);
    const auto *nsg0_980 = buffer.data(nsg0 + 980);
    const auto *nsg0_985 = buffer.data(nsg0 + 985);
    const auto *nsg0_986 = buffer.data(nsg0 + 986);
    const auto *nsg0_987 = buffer.data(nsg0 + 987);
    const auto *nsg0_989 = buffer.data(nsg0 + 989);

    const auto *nsf_530 = buffer.data(nsf + 530);
    const auto *nsf_536 = buffer.data(nsf + 536);
    const auto *nsf_540 = buffer.data(nsf + 540);
    const auto *nsf_542 = buffer.data(nsf + 542);
    const auto *nsf_549 = buffer.data(nsf + 549);
    const auto *nsf_556 = buffer.data(nsf + 556);
    const auto *nsf_557 = buffer.data(nsf + 557);
    const auto *nsf_559 = buffer.data(nsf + 559);
    const auto *nsf_566 = buffer.data(nsf + 566);
    const auto *nsf_569 = buffer.data(nsf + 569);
    const auto *nsf_576 = buffer.data(nsf + 576);
    const auto *nsf_578 = buffer.data(nsf + 578);
    const auto *nsf_579 = buffer.data(nsf + 579);
    const auto *nsf_586 = buffer.data(nsf + 586);
    const auto *nsf_588 = buffer.data(nsf + 588);
    const auto *nsf_589 = buffer.data(nsf + 589);
    const auto *nsf_596 = buffer.data(nsf + 596);
    const auto *nsf_598 = buffer.data(nsf + 598);
    const auto *nsf_599 = buffer.data(nsf + 599);
    const auto *nsf_606 = buffer.data(nsf + 606);
    const auto *nsf_608 = buffer.data(nsf + 608);
    const auto *nsf_609 = buffer.data(nsf + 609);
    const auto *nsf_643 = buffer.data(nsf + 643);
    const auto *nsf_646 = buffer.data(nsf + 646);
    const auto *nsf_647 = buffer.data(nsf + 647);
    const auto *nsf_648 = buffer.data(nsf + 648);
    const auto *nsf_649 = buffer.data(nsf + 649);
    const auto *nsf_650 = buffer.data(nsf + 650);
    const auto *nsf_655 = buffer.data(nsf + 655);
    const auto *nsf_656 = buffer.data(nsf + 656);
    const auto *nsf_657 = buffer.data(nsf + 657);
    const auto *nsf_659 = buffer.data(nsf + 659);

    const auto *nsg1_815 = buffer.data(nsg1 + 815);
    const auto *nsg1_825 = buffer.data(nsg1 + 825);
    const auto *nsg1_826 = buffer.data(nsg1 + 826);
    const auto *nsg1_828 = buffer.data(nsg1 + 828);
    const auto *nsg1_835 = buffer.data(nsg1 + 835);
    const auto *nsg1_837 = buffer.data(nsg1 + 837);
    const auto *nsg1_963 = buffer.data(nsg1 + 963);
    const auto *nsg1_970 = buffer.data(nsg1 + 970);
    const auto *nsg1_972 = buffer.data(nsg1 + 972);
    const auto *nsg1_974 = buffer.data(nsg1 + 974);
    const auto *nsg1_975 = buffer.data(nsg1 + 975);
    const auto *nsg1_980 = buffer.data(nsg1 + 980);
    const auto *nsg1_985 = buffer.data(nsg1 + 985);
    const auto *nsg1_986 = buffer.data(nsg1 + 986);
    const auto *nsg1_987 = buffer.data(nsg1 + 987);
    const auto *nsg1_989 = buffer.data(nsg1 + 989);

    const auto *osd0_390 = buffer.data(osd0 + 390);
    const auto *osd0_396 = buffer.data(osd0 + 396);
    const auto *osd0_397 = buffer.data(osd0 + 397);
    const auto *osd0_399 = buffer.data(osd0 + 399);
    const auto *osd0_401 = buffer.data(osd0 + 401);
    const auto *osd0_404 = buffer.data(osd0 + 404);
    const auto *osd0_406 = buffer.data(osd0 + 406);
    const auto *osd0_407 = buffer.data(osd0 + 407);
    const auto *osd0_408 = buffer.data(osd0 + 408);
    const auto *osd0_409 = buffer.data(osd0 + 409);
    const auto *osd0_410 = buffer.data(osd0 + 410);
    const auto *osd0_411 = buffer.data(osd0 + 411);
    const auto *osd0_412 = buffer.data(osd0 + 412);
    const auto *osd0_413 = buffer.data(osd0 + 413);
    const auto *osd0_414 = buffer.data(osd0 + 414);
    const auto *osd0_415 = buffer.data(osd0 + 415);
    const auto *osd0_416 = buffer.data(osd0 + 416);
    const auto *osd0_417 = buffer.data(osd0 + 417);
    const auto *osd0_418 = buffer.data(osd0 + 418);
    const auto *osd0_419 = buffer.data(osd0 + 419);
    const auto *osd0_420 = buffer.data(osd0 + 420);
    const auto *osd0_421 = buffer.data(osd0 + 421);
    const auto *osd0_422 = buffer.data(osd0 + 422);
    const auto *osd0_423 = buffer.data(osd0 + 423);
    const auto *osd0_424 = buffer.data(osd0 + 424);
    const auto *osd0_425 = buffer.data(osd0 + 425);
    const auto *osd0_426 = buffer.data(osd0 + 426);
    const auto *osd0_427 = buffer.data(osd0 + 427);
    const auto *osd0_428 = buffer.data(osd0 + 428);
    const auto *osd0_429 = buffer.data(osd0 + 429);
    const auto *osd0_430 = buffer.data(osd0 + 430);
    const auto *osd0_431 = buffer.data(osd0 + 431);
    const auto *osd0_432 = buffer.data(osd0 + 432);
    const auto *osd0_433 = buffer.data(osd0 + 433);
    const auto *osd0_434 = buffer.data(osd0 + 434);
    const auto *osd0_435 = buffer.data(osd0 + 435);
    const auto *osd0_436 = buffer.data(osd0 + 436);
    const auto *osd0_437 = buffer.data(osd0 + 437);

    const auto *osd1_390 = buffer.data(osd1 + 390);
    const auto *osd1_396 = buffer.data(osd1 + 396);
    const auto *osd1_397 = buffer.data(osd1 + 397);
    const auto *osd1_399 = buffer.data(osd1 + 399);
    const auto *osd1_401 = buffer.data(osd1 + 401);
    const auto *osd1_404 = buffer.data(osd1 + 404);
    const auto *osd1_406 = buffer.data(osd1 + 406);
    const auto *osd1_407 = buffer.data(osd1 + 407);
    const auto *osd1_408 = buffer.data(osd1 + 408);
    const auto *osd1_409 = buffer.data(osd1 + 409);
    const auto *osd1_410 = buffer.data(osd1 + 410);
    const auto *osd1_411 = buffer.data(osd1 + 411);
    const auto *osd1_412 = buffer.data(osd1 + 412);
    const auto *osd1_413 = buffer.data(osd1 + 413);
    const auto *osd1_414 = buffer.data(osd1 + 414);
    const auto *osd1_415 = buffer.data(osd1 + 415);
    const auto *osd1_416 = buffer.data(osd1 + 416);
    const auto *osd1_417 = buffer.data(osd1 + 417);
    const auto *osd1_418 = buffer.data(osd1 + 418);
    const auto *osd1_419 = buffer.data(osd1 + 419);
    const auto *osd1_420 = buffer.data(osd1 + 420);
    const auto *osd1_421 = buffer.data(osd1 + 421);
    const auto *osd1_422 = buffer.data(osd1 + 422);
    const auto *osd1_423 = buffer.data(osd1 + 423);
    const auto *osd1_424 = buffer.data(osd1 + 424);
    const auto *osd1_425 = buffer.data(osd1 + 425);
    const auto *osd1_426 = buffer.data(osd1 + 426);
    const auto *osd1_427 = buffer.data(osd1 + 427);
    const auto *osd1_428 = buffer.data(osd1 + 428);
    const auto *osd1_429 = buffer.data(osd1 + 429);
    const auto *osd1_430 = buffer.data(osd1 + 430);
    const auto *osd1_431 = buffer.data(osd1 + 431);
    const auto *osd1_432 = buffer.data(osd1 + 432);
    const auto *osd1_433 = buffer.data(osd1 + 433);
    const auto *osd1_434 = buffer.data(osd1 + 434);
    const auto *osd1_435 = buffer.data(osd1 + 435);
    const auto *osd1_436 = buffer.data(osd1 + 436);
    const auto *osd1_437 = buffer.data(osd1 + 437);

    const auto *osf_640 = buffer.data(osf + 640);
    const auto *osf_642 = buffer.data(osf + 642);
    const auto *osf_646 = buffer.data(osf + 646);
    const auto *osf_647 = buffer.data(osf + 647);
    const auto *osf_648 = buffer.data(osf + 648);
    const auto *osf_649 = buffer.data(osf + 649);
    const auto *osf_650 = buffer.data(osf + 650);
    const auto *osf_651 = buffer.data(osf + 651);
    const auto *osf_652 = buffer.data(osf + 652);
    const auto *osf_655 = buffer.data(osf + 655);
    const auto *osf_656 = buffer.data(osf + 656);
    const auto *osf_657 = buffer.data(osf + 657);
    const auto *osf_659 = buffer.data(osf + 659);
    const auto *osf_660 = buffer.data(osf + 660);
    const auto *osf_661 = buffer.data(osf + 661);
    const auto *osf_663 = buffer.data(osf + 663);
    const auto *osf_665 = buffer.data(osf + 665);
    const auto *osf_666 = buffer.data(osf + 666);
    const auto *osf_667 = buffer.data(osf + 667);
    const auto *osf_668 = buffer.data(osf + 668);
    const auto *osf_669 = buffer.data(osf + 669);
    const auto *osf_672 = buffer.data(osf + 672);
    const auto *osf_674 = buffer.data(osf + 674);
    const auto *osf_675 = buffer.data(osf + 675);
    const auto *osf_676 = buffer.data(osf + 676);
    const auto *osf_677 = buffer.data(osf + 677);
    const auto *osf_678 = buffer.data(osf + 678);
    const auto *osf_679 = buffer.data(osf + 679);
    const auto *osf_680 = buffer.data(osf + 680);
    const auto *osf_681 = buffer.data(osf + 681);
    const auto *osf_682 = buffer.data(osf + 682);
    const auto *osf_683 = buffer.data(osf + 683);
    const auto *osf_684 = buffer.data(osf + 684);
    const auto *osf_685 = buffer.data(osf + 685);
    const auto *osf_686 = buffer.data(osf + 686);
    const auto *osf_687 = buffer.data(osf + 687);
    const auto *osf_688 = buffer.data(osf + 688);
    const auto *osf_689 = buffer.data(osf + 689);
    const auto *osf_690 = buffer.data(osf + 690);
    const auto *osf_691 = buffer.data(osf + 691);
    const auto *osf_692 = buffer.data(osf + 692);
    const auto *osf_693 = buffer.data(osf + 693);
    const auto *osf_694 = buffer.data(osf + 694);
    const auto *osf_695 = buffer.data(osf + 695);
    const auto *osf_696 = buffer.data(osf + 696);
    const auto *osf_697 = buffer.data(osf + 697);
    const auto *osf_698 = buffer.data(osf + 698);
    const auto *osf_699 = buffer.data(osf + 699);
    const auto *osf_700 = buffer.data(osf + 700);
    const auto *osf_701 = buffer.data(osf + 701);
    const auto *osf_702 = buffer.data(osf + 702);
    const auto *osf_703 = buffer.data(osf + 703);
    const auto *osf_704 = buffer.data(osf + 704);
    const auto *osf_705 = buffer.data(osf + 705);
    const auto *osf_706 = buffer.data(osf + 706);
    const auto *osf_707 = buffer.data(osf + 707);
    const auto *osf_708 = buffer.data(osf + 708);
    const auto *osf_709 = buffer.data(osf + 709);
    const auto *osf_710 = buffer.data(osf + 710);
    const auto *osf_711 = buffer.data(osf + 711);
    const auto *osf_712 = buffer.data(osf + 712);
    const auto *osf_713 = buffer.data(osf + 713);
    const auto *osf_714 = buffer.data(osf + 714);
    const auto *osf_715 = buffer.data(osf + 715);
    const auto *osf_716 = buffer.data(osf + 716);
    const auto *osf_717 = buffer.data(osf + 717);
    const auto *osf_718 = buffer.data(osf + 718);
    const auto *osf_719 = buffer.data(osf + 719);
    const auto *osf_720 = buffer.data(osf + 720);
    const auto *osf_721 = buffer.data(osf + 721);
    const auto *osf_722 = buffer.data(osf + 722);
    const auto *osf_723 = buffer.data(osf + 723);
    const auto *osf_724 = buffer.data(osf + 724);
    const auto *osf_725 = buffer.data(osf + 725);
    const auto *osf_726 = buffer.data(osf + 726);

#pragma omp simd aligned(t_962, t_963, t_964, pa_x, pc_x, pc_y, pc_z, nsg0_963, nsf_530, \
                         nsf_542, nsf_643, nsg1_963, osf_640, osf_642 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_962[k] = f_12 * nsf_530[k]
                   + f_3 * pc_z[k] * osf_640[k];

        t_963[k] = pa_x[k] * nsg0_963[k]
                   + f_8 * nsf_643[k]
                   - f_6 * pc_x[k] * nsg1_963[k];

        t_964[k] = f_7 * nsf_542[k]
                   + f_3 * pc_y[k] * osf_642[k];
    }

#pragma omp simd aligned(t_965, t_966, t_967, t_968, pa_y, pc_x, pc_y, nsg0_815, nsf_646, \
                         nsf_647, nsf_648, nsg1_815, osf_646, osf_647, \
                         osf_648 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_965[k] = pa_y[k] * nsg0_815[k]
                   - f_6 * pc_y[k] * nsg1_815[k];

        t_966[k] = f_7 * nsf_646[k]
                   + f_3 * pc_x[k] * osf_646[k];

        t_967[k] = f_7 * nsf_647[k]
                   + f_3 * pc_x[k] * osf_647[k];

        t_968[k] = f_7 * nsf_648[k]
                   + f_3 * pc_x[k] * osf_648[k];
    }

#pragma omp simd aligned(t_969, t_970, t_971, t_972, pa_x, pc_x, pc_z, nsg0_970, nsg0_972, \
                         nsf_536, nsf_649, nsg1_970, nsg1_972, osf_646, \
                         osf_649 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_969[k] = f_7 * nsf_649[k]
                   + f_3 * pc_x[k] * osf_649[k];

        t_970[k] = pa_x[k] * nsg0_970[k]
                   - f_6 * pc_x[k] * nsg1_970[k];

        t_971[k] = f_12 * nsf_536[k]
                   + f_3 * pc_z[k] * osf_646[k];

        t_972[k] = pa_x[k] * nsg0_972[k]
                   - f_6 * pc_x[k] * nsg1_972[k];
    }

#pragma omp simd aligned(t_973, t_974, t_975, t_976, pa_x, pc_x, pc_y, nsg0_974, nsg0_975, \
                         nsf_549, nsf_650, nsg1_974, nsg1_975, osf_649, \
                         osf_650 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_973[k] = f_7 * nsf_549[k]
                   + f_3 * pc_y[k] * osf_649[k];

        t_974[k] = pa_x[k] * nsg0_974[k]
                   - f_6 * pc_x[k] * nsg1_974[k];

        t_975[k] = pa_x[k] * nsg0_975[k]
                   + f_16 * nsf_650[k]
                   - f_6 * pc_x[k] * nsg1_975[k];

        t_976[k] = f_3 * pc_y[k] * osf_650[k];
    }

#pragma omp simd aligned(t_977, t_978, t_979, pc_y, pc_z, nsf_540, osd0_390, osd1_390, \
                         osf_650, osf_651, osf_652 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_977[k] = f_9 * nsf_540[k]
                   + f_3 * pc_z[k] * osf_650[k];

        t_978[k] = f_4 * osd0_390[k]
                   - f_5 * osd1_390[k]
                   + f_3 * pc_y[k] * osf_651[k];

        t_979[k] = f_3 * pc_y[k] * osf_652[k];
    }

#pragma omp simd aligned(t_980, t_981, t_982, t_983, pa_x, pc_x, pc_y, nsg0_980, nsf_655, \
                         nsf_656, nsf_657, nsg1_980, osf_655, osf_656, \
                         osf_657 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_980[k] = pa_x[k] * nsg0_980[k]
                   + f_8 * nsf_655[k]
                   - f_6 * pc_x[k] * nsg1_980[k];

        t_981[k] = f_7 * nsf_656[k]
                   + f_3 * pc_x[k] * osf_656[k];

        t_982[k] = f_7 * nsf_657[k]
                   + f_3 * pc_x[k] * osf_657[k];

        t_983[k] = f_3 * pc_y[k] * osf_655[k];
    }

#pragma omp simd aligned(t_984, t_985, t_986, t_987, t_988, pa_x, pc_x, pc_y, nsg0_985, \
                         nsg0_986, nsg0_987, nsf_659, nsg1_985, nsg1_986, nsg1_987, \
                         osf_659 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_984[k] = f_7 * nsf_659[k]
                   + f_3 * pc_x[k] * osf_659[k];

        t_985[k] = pa_x[k] * nsg0_985[k]
                   - f_6 * pc_x[k] * nsg1_985[k];

        t_986[k] = pa_x[k] * nsg0_986[k]
                   - f_6 * pc_x[k] * nsg1_986[k];

        t_987[k] = pa_x[k] * nsg0_987[k]
                   - f_6 * pc_x[k] * nsg1_987[k];

        t_988[k] = f_3 * pc_y[k] * osf_659[k];
    }

#pragma omp simd aligned(t_989, t_990, t_991, t_992, pa_x, pc_x, pc_z, nsg0_989, nsg1_989, \
                         osd0_396, osd0_397, osd1_396, osd1_397, osf_660, \
                         osf_661 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_989[k] = pa_x[k] * nsg0_989[k]
                   - f_6 * pc_x[k] * nsg1_989[k];

        t_990[k] = f_1 * osd0_396[k]
                   - f_2 * osd1_396[k]
                   + f_3 * pc_x[k] * osf_660[k];

        t_991[k] = f_10 * osd0_397[k]
                   - f_11 * osd1_397[k]
                   + f_3 * pc_x[k] * osf_661[k];

        t_992[k] = f_3 * pc_z[k] * osf_660[k];
    }

#pragma omp simd aligned(t_993, t_994, t_995, t_996, t_997, pc_x, pc_z, osd0_399, osd0_401, \
                         osd1_399, osd1_401, osf_661, osf_663, osf_665, osf_666, \
                         osf_667 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_993[k] = f_4 * osd0_399[k]
                   - f_5 * osd1_399[k]
                   + f_3 * pc_x[k] * osf_663[k];

        t_994[k] = f_3 * pc_z[k] * osf_661[k];

        t_995[k] = f_4 * osd0_401[k]
                   - f_5 * osd1_401[k]
                   + f_3 * pc_x[k] * osf_665[k];

        t_996[k] = f_3 * pc_x[k] * osf_666[k];

        t_997[k] = f_3 * pc_x[k] * osf_667[k];
    }

#pragma omp simd aligned(t_998, t_999, t_1000, t_1001, t_1002, pc_x, pc_y, pc_z, nsf_556, \
                         osd0_399, osd1_399, osf_666, osf_667, osf_668, \
                         osf_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_998[k] = f_3 * pc_x[k] * osf_668[k];

        t_999[k] = f_3 * pc_x[k] * osf_669[k];

        t_1000[k] = f_0 * nsf_556[k]
                    + f_1 * osd0_399[k]
                    - f_2 * osd1_399[k]
                    + f_3 * pc_y[k] * osf_666[k];

        t_1001[k] = f_3 * pc_z[k] * osf_666[k];

        t_1002[k] = f_4 * osd0_399[k]
                    - f_5 * osd1_399[k]
                    + f_3 * pc_z[k] * osf_667[k];
    }

#pragma omp simd aligned(t_1003, t_1004, t_1005, t_1006, pa_z, pc_y, pc_z, nsg0_825, nsg0_826, \
                         nsf_559, nsg1_825, nsg1_826, osd0_401, osd1_401, \
                         osf_669 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1003[k] = f_0 * nsf_559[k]
                    + f_3 * pc_y[k] * osf_669[k];

        t_1004[k] = f_1 * osd0_401[k]
                    - f_2 * osd1_401[k]
                    + f_3 * pc_z[k] * osf_669[k];

        t_1005[k] = pa_z[k] * nsg0_825[k]
                    - f_6 * pc_z[k] * nsg1_825[k];

        t_1006[k] = pa_z[k] * nsg0_826[k]
                    - f_6 * pc_z[k] * nsg1_826[k];
    }

#pragma omp simd aligned(t_1007, t_1008, t_1009, pa_z, pc_x, pc_z, nsg0_828, nsg1_828, \
                         osd0_404, osd0_406, osd1_404, osd1_406, osf_672, \
                         osf_674 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1007[k] = f_10 * osd0_404[k]
                    - f_11 * osd1_404[k]
                    + f_3 * pc_x[k] * osf_672[k];

        t_1008[k] = pa_z[k] * nsg0_828[k]
                    - f_6 * pc_z[k] * nsg1_828[k];

        t_1009[k] = f_4 * osd0_406[k]
                    - f_5 * osd1_406[k]
                    + f_3 * pc_x[k] * osf_674[k];
    }

#pragma omp simd aligned(t_1010, t_1011, t_1012, t_1013, t_1014, pc_x, osd0_407, osd1_407, \
                         osf_675, osf_676, osf_677, osf_678, osf_679 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1010[k] = f_4 * osd0_407[k]
                    - f_5 * osd1_407[k]
                    + f_3 * pc_x[k] * osf_675[k];

        t_1011[k] = f_3 * pc_x[k] * osf_676[k];

        t_1012[k] = f_3 * pc_x[k] * osf_677[k];

        t_1013[k] = f_3 * pc_x[k] * osf_678[k];

        t_1014[k] = f_3 * pc_x[k] * osf_679[k];
    }

#pragma omp simd aligned(t_1015, t_1016, t_1017, t_1018, pa_z, pc_y, pc_z, nsg0_835, nsg0_837, \
                         nsf_556, nsf_557, nsf_569, nsg1_835, nsg1_837, osf_676, \
                         osf_679 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1015[k] = pa_z[k] * nsg0_835[k]
                    - f_6 * pc_z[k] * nsg1_835[k];

        t_1016[k] = f_7 * nsf_556[k]
                    + f_3 * pc_z[k] * osf_676[k];

        t_1017[k] = pa_z[k] * nsg0_837[k]
                    + f_8 * nsf_557[k]
                    - f_6 * pc_z[k] * nsg1_837[k];

        t_1018[k] = f_9 * nsf_569[k]
                    + f_3 * pc_y[k] * osf_679[k];
    }

#pragma omp simd aligned(t_1019, t_1020, t_1021, pc_x, pc_z, nsf_559, osd0_407, osd0_408, \
                         osd0_409, osd1_407, osd1_408, osd1_409, osf_679, osf_680, \
                         osf_681 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1019[k] = f_7 * nsf_559[k]
                    + f_1 * osd0_407[k]
                    - f_2 * osd1_407[k]
                    + f_3 * pc_z[k] * osf_679[k];

        t_1020[k] = f_1 * osd0_408[k]
                    - f_2 * osd1_408[k]
                    + f_3 * pc_x[k] * osf_680[k];

        t_1021[k] = f_10 * osd0_409[k]
                    - f_11 * osd1_409[k]
                    + f_3 * pc_x[k] * osf_681[k];
    }

#pragma omp simd aligned(t_1022, t_1023, t_1024, pc_x, osd0_410, osd0_411, osd0_412, osd1_410, \
                         osd1_411, osd1_412, osf_682, osf_683, \
                         osf_684 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1022[k] = f_10 * osd0_410[k]
                    - f_11 * osd1_410[k]
                    + f_3 * pc_x[k] * osf_682[k];

        t_1023[k] = f_4 * osd0_411[k]
                    - f_5 * osd1_411[k]
                    + f_3 * pc_x[k] * osf_683[k];

        t_1024[k] = f_4 * osd0_412[k]
                    - f_5 * osd1_412[k]
                    + f_3 * pc_x[k] * osf_684[k];
    }

#pragma omp simd aligned(t_1025, t_1026, t_1027, t_1028, t_1029, pc_x, osd0_413, osd1_413, \
                         osf_685, osf_686, osf_687, osf_688, osf_689 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1025[k] = f_4 * osd0_413[k]
                    - f_5 * osd1_413[k]
                    + f_3 * pc_x[k] * osf_685[k];

        t_1026[k] = f_3 * pc_x[k] * osf_686[k];

        t_1027[k] = f_3 * pc_x[k] * osf_687[k];

        t_1028[k] = f_3 * pc_x[k] * osf_688[k];

        t_1029[k] = f_3 * pc_x[k] * osf_689[k];
    }

#pragma omp simd aligned(t_1030, t_1031, t_1032, pc_y, pc_z, nsf_566, nsf_576, nsf_578, \
                         osd0_411, osd0_413, osd1_411, osd1_413, osf_686, \
                         osf_688 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1030[k] = f_12 * nsf_576[k]
                    + f_1 * osd0_411[k]
                    - f_2 * osd1_411[k]
                    + f_3 * pc_y[k] * osf_686[k];

        t_1031[k] = f_8 * nsf_566[k]
                    + f_3 * pc_z[k] * osf_686[k];

        t_1032[k] = f_12 * nsf_578[k]
                    + f_4 * osd0_413[k]
                    - f_5 * osd1_413[k]
                    + f_3 * pc_y[k] * osf_688[k];
    }

#pragma omp simd aligned(t_1033, t_1034, t_1035, pc_x, pc_y, pc_z, nsf_569, nsf_579, osd0_413, \
                         osd0_414, osd1_413, osd1_414, osf_689, \
                         osf_690 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1033[k] = f_12 * nsf_579[k]
                    + f_3 * pc_y[k] * osf_689[k];

        t_1034[k] = f_8 * nsf_569[k]
                    + f_1 * osd0_413[k]
                    - f_2 * osd1_413[k]
                    + f_3 * pc_z[k] * osf_689[k];

        t_1035[k] = f_1 * osd0_414[k]
                    - f_2 * osd1_414[k]
                    + f_3 * pc_x[k] * osf_690[k];
    }

#pragma omp simd aligned(t_1036, t_1037, t_1038, pc_x, osd0_415, osd0_416, osd0_417, osd1_415, \
                         osd1_416, osd1_417, osf_691, osf_692, \
                         osf_693 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1036[k] = f_10 * osd0_415[k]
                    - f_11 * osd1_415[k]
                    + f_3 * pc_x[k] * osf_691[k];

        t_1037[k] = f_10 * osd0_416[k]
                    - f_11 * osd1_416[k]
                    + f_3 * pc_x[k] * osf_692[k];

        t_1038[k] = f_4 * osd0_417[k]
                    - f_5 * osd1_417[k]
                    + f_3 * pc_x[k] * osf_693[k];
    }

#pragma omp simd aligned(t_1039, t_1040, t_1041, t_1042, t_1043, pc_x, osd0_418, osd0_419, \
                         osd1_418, osd1_419, osf_694, osf_695, osf_696, osf_697, \
                         osf_698 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1039[k] = f_4 * osd0_418[k]
                    - f_5 * osd1_418[k]
                    + f_3 * pc_x[k] * osf_694[k];

        t_1040[k] = f_4 * osd0_419[k]
                    - f_5 * osd1_419[k]
                    + f_3 * pc_x[k] * osf_695[k];

        t_1041[k] = f_3 * pc_x[k] * osf_696[k];

        t_1042[k] = f_3 * pc_x[k] * osf_697[k];

        t_1043[k] = f_3 * pc_x[k] * osf_698[k];
    }

#pragma omp simd aligned(t_1044, t_1045, t_1046, pc_x, pc_y, pc_z, nsf_576, nsf_586, osd0_417, \
                         osd1_417, osf_696, osf_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1044[k] = f_3 * pc_x[k] * osf_699[k];

        t_1045[k] = f_13 * nsf_586[k]
                    + f_1 * osd0_417[k]
                    - f_2 * osd1_417[k]
                    + f_3 * pc_y[k] * osf_696[k];

        t_1046[k] = f_14 * nsf_576[k]
                    + f_3 * pc_z[k] * osf_696[k];
    }

#pragma omp simd aligned(t_1047, t_1048, t_1049, pc_y, pc_z, nsf_579, nsf_588, nsf_589, \
                         osd0_419, osd1_419, osf_698, osf_699 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1047[k] = f_13 * nsf_588[k]
                    + f_4 * osd0_419[k]
                    - f_5 * osd1_419[k]
                    + f_3 * pc_y[k] * osf_698[k];

        t_1048[k] = f_13 * nsf_589[k]
                    + f_3 * pc_y[k] * osf_699[k];

        t_1049[k] = f_14 * nsf_579[k]
                    + f_1 * osd0_419[k]
                    - f_2 * osd1_419[k]
                    + f_3 * pc_z[k] * osf_699[k];
    }

#pragma omp simd aligned(t_1050, t_1051, t_1052, pc_x, osd0_420, osd0_421, osd0_422, osd1_420, \
                         osd1_421, osd1_422, osf_700, osf_701, \
                         osf_702 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1050[k] = f_1 * osd0_420[k]
                    - f_2 * osd1_420[k]
                    + f_3 * pc_x[k] * osf_700[k];

        t_1051[k] = f_10 * osd0_421[k]
                    - f_11 * osd1_421[k]
                    + f_3 * pc_x[k] * osf_701[k];

        t_1052[k] = f_10 * osd0_422[k]
                    - f_11 * osd1_422[k]
                    + f_3 * pc_x[k] * osf_702[k];
    }

#pragma omp simd aligned(t_1053, t_1054, t_1055, t_1056, pc_x, osd0_423, osd0_424, osd0_425, \
                         osd1_423, osd1_424, osd1_425, osf_703, osf_704, osf_705, \
                         osf_706 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1053[k] = f_4 * osd0_423[k]
                    - f_5 * osd1_423[k]
                    + f_3 * pc_x[k] * osf_703[k];

        t_1054[k] = f_4 * osd0_424[k]
                    - f_5 * osd1_424[k]
                    + f_3 * pc_x[k] * osf_704[k];

        t_1055[k] = f_4 * osd0_425[k]
                    - f_5 * osd1_425[k]
                    + f_3 * pc_x[k] * osf_705[k];

        t_1056[k] = f_3 * pc_x[k] * osf_706[k];
    }

#pragma omp simd aligned(t_1057, t_1058, t_1059, t_1060, t_1061, pc_x, pc_y, pc_z, nsf_586, \
                         nsf_596, osd0_423, osd1_423, osf_706, osf_707, osf_708, \
                         osf_709 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1057[k] = f_3 * pc_x[k] * osf_707[k];

        t_1058[k] = f_3 * pc_x[k] * osf_708[k];

        t_1059[k] = f_3 * pc_x[k] * osf_709[k];

        t_1060[k] = f_15 * nsf_596[k]
                    + f_1 * osd0_423[k]
                    - f_2 * osd1_423[k]
                    + f_3 * pc_y[k] * osf_706[k];

        t_1061[k] = f_16 * nsf_586[k]
                    + f_3 * pc_z[k] * osf_706[k];
    }

#pragma omp simd aligned(t_1062, t_1063, t_1064, pc_y, pc_z, nsf_589, nsf_598, nsf_599, \
                         osd0_425, osd1_425, osf_708, osf_709 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1062[k] = f_15 * nsf_598[k]
                    + f_4 * osd0_425[k]
                    - f_5 * osd1_425[k]
                    + f_3 * pc_y[k] * osf_708[k];

        t_1063[k] = f_15 * nsf_599[k]
                    + f_3 * pc_y[k] * osf_709[k];

        t_1064[k] = f_16 * nsf_589[k]
                    + f_1 * osd0_425[k]
                    - f_2 * osd1_425[k]
                    + f_3 * pc_z[k] * osf_709[k];
    }

#pragma omp simd aligned(t_1065, t_1066, t_1067, pc_x, osd0_426, osd0_427, osd0_428, osd1_426, \
                         osd1_427, osd1_428, osf_710, osf_711, \
                         osf_712 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1065[k] = f_1 * osd0_426[k]
                    - f_2 * osd1_426[k]
                    + f_3 * pc_x[k] * osf_710[k];

        t_1066[k] = f_10 * osd0_427[k]
                    - f_11 * osd1_427[k]
                    + f_3 * pc_x[k] * osf_711[k];

        t_1067[k] = f_10 * osd0_428[k]
                    - f_11 * osd1_428[k]
                    + f_3 * pc_x[k] * osf_712[k];
    }

#pragma omp simd aligned(t_1068, t_1069, t_1070, t_1071, pc_x, osd0_429, osd0_430, osd0_431, \
                         osd1_429, osd1_430, osd1_431, osf_713, osf_714, osf_715, \
                         osf_716 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1068[k] = f_4 * osd0_429[k]
                    - f_5 * osd1_429[k]
                    + f_3 * pc_x[k] * osf_713[k];

        t_1069[k] = f_4 * osd0_430[k]
                    - f_5 * osd1_430[k]
                    + f_3 * pc_x[k] * osf_714[k];

        t_1070[k] = f_4 * osd0_431[k]
                    - f_5 * osd1_431[k]
                    + f_3 * pc_x[k] * osf_715[k];

        t_1071[k] = f_3 * pc_x[k] * osf_716[k];
    }

#pragma omp simd aligned(t_1072, t_1073, t_1074, t_1075, t_1076, pc_x, pc_y, pc_z, nsf_596, \
                         nsf_606, osd0_429, osd1_429, osf_716, osf_717, osf_718, \
                         osf_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1072[k] = f_3 * pc_x[k] * osf_717[k];

        t_1073[k] = f_3 * pc_x[k] * osf_718[k];

        t_1074[k] = f_3 * pc_x[k] * osf_719[k];

        t_1075[k] = f_17 * nsf_606[k]
                    + f_1 * osd0_429[k]
                    - f_2 * osd1_429[k]
                    + f_3 * pc_y[k] * osf_716[k];

        t_1076[k] = f_18 * nsf_596[k]
                    + f_3 * pc_z[k] * osf_716[k];
    }

#pragma omp simd aligned(t_1077, t_1078, t_1079, pc_y, pc_z, nsf_599, nsf_608, nsf_609, \
                         osd0_431, osd1_431, osf_718, osf_719 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1077[k] = f_17 * nsf_608[k]
                    + f_4 * osd0_431[k]
                    - f_5 * osd1_431[k]
                    + f_3 * pc_y[k] * osf_718[k];

        t_1078[k] = f_17 * nsf_609[k]
                    + f_3 * pc_y[k] * osf_719[k];

        t_1079[k] = f_18 * nsf_599[k]
                    + f_1 * osd0_431[k]
                    - f_2 * osd1_431[k]
                    + f_3 * pc_z[k] * osf_719[k];
    }

#pragma omp simd aligned(t_1080, t_1081, t_1082, pc_x, osd0_432, osd0_433, osd0_434, osd1_432, \
                         osd1_433, osd1_434, osf_720, osf_721, \
                         osf_722 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1080[k] = f_1 * osd0_432[k]
                    - f_2 * osd1_432[k]
                    + f_3 * pc_x[k] * osf_720[k];

        t_1081[k] = f_10 * osd0_433[k]
                    - f_11 * osd1_433[k]
                    + f_3 * pc_x[k] * osf_721[k];

        t_1082[k] = f_10 * osd0_434[k]
                    - f_11 * osd1_434[k]
                    + f_3 * pc_x[k] * osf_722[k];
    }

#pragma omp simd aligned(t_1083, t_1084, t_1085, t_1086, pc_x, osd0_435, osd0_436, osd0_437, \
                         osd1_435, osd1_436, osd1_437, osf_723, osf_724, osf_725, \
                         osf_726 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1083[k] = f_4 * osd0_435[k]
                    - f_5 * osd1_435[k]
                    + f_3 * pc_x[k] * osf_723[k];

        t_1084[k] = f_4 * osd0_436[k]
                    - f_5 * osd1_436[k]
                    + f_3 * pc_x[k] * osf_724[k];

        t_1085[k] = f_4 * osd0_437[k]
                    - f_5 * osd1_437[k]
                    + f_3 * pc_x[k] * osf_725[k];

        t_1086[k] = f_3 * pc_x[k] * osf_726[k];
    }
}

static auto
compute_prim_osg_three_center_electron_repulsion_0_piece9(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nsg0,
                                                          const size_t nsf, const size_t nsg1,
                                                          const size_t osd0, const size_t osd1,
                                                          const size_t osf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 5.0 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 4.5 / q;
    const auto f_13 = 4.0 / q;
    const auto f_14 = 1.5 / q;
    const auto f_15 = 3.5 / q;
    const auto f_16 = 2.0 / q;
    const auto f_17 = 3.0 / q;
    const auto f_18 = 2.5 / q;

    auto *t_1087 = buffer.data(target + 1087);
    auto *t_1088 = buffer.data(target + 1088);
    auto *t_1089 = buffer.data(target + 1089);
    auto *t_1090 = buffer.data(target + 1090);
    auto *t_1091 = buffer.data(target + 1091);
    auto *t_1092 = buffer.data(target + 1092);
    auto *t_1093 = buffer.data(target + 1093);
    auto *t_1094 = buffer.data(target + 1094);
    auto *t_1095 = buffer.data(target + 1095);
    auto *t_1096 = buffer.data(target + 1096);
    auto *t_1097 = buffer.data(target + 1097);
    auto *t_1098 = buffer.data(target + 1098);
    auto *t_1099 = buffer.data(target + 1099);
    auto *t_1100 = buffer.data(target + 1100);
    auto *t_1101 = buffer.data(target + 1101);
    auto *t_1102 = buffer.data(target + 1102);
    auto *t_1103 = buffer.data(target + 1103);
    auto *t_1104 = buffer.data(target + 1104);
    auto *t_1105 = buffer.data(target + 1105);
    auto *t_1106 = buffer.data(target + 1106);
    auto *t_1107 = buffer.data(target + 1107);
    auto *t_1108 = buffer.data(target + 1108);
    auto *t_1109 = buffer.data(target + 1109);
    auto *t_1110 = buffer.data(target + 1110);
    auto *t_1111 = buffer.data(target + 1111);
    auto *t_1112 = buffer.data(target + 1112);
    auto *t_1113 = buffer.data(target + 1113);
    auto *t_1114 = buffer.data(target + 1114);
    auto *t_1115 = buffer.data(target + 1115);
    auto *t_1116 = buffer.data(target + 1116);
    auto *t_1117 = buffer.data(target + 1117);
    auto *t_1118 = buffer.data(target + 1118);
    auto *t_1119 = buffer.data(target + 1119);
    auto *t_1120 = buffer.data(target + 1120);
    auto *t_1121 = buffer.data(target + 1121);
    auto *t_1122 = buffer.data(target + 1122);
    auto *t_1123 = buffer.data(target + 1123);
    auto *t_1124 = buffer.data(target + 1124);
    auto *t_1125 = buffer.data(target + 1125);
    auto *t_1126 = buffer.data(target + 1126);
    auto *t_1127 = buffer.data(target + 1127);
    auto *t_1128 = buffer.data(target + 1128);
    auto *t_1129 = buffer.data(target + 1129);
    auto *t_1130 = buffer.data(target + 1130);
    auto *t_1131 = buffer.data(target + 1131);
    auto *t_1132 = buffer.data(target + 1132);
    auto *t_1133 = buffer.data(target + 1133);
    auto *t_1134 = buffer.data(target + 1134);
    auto *t_1135 = buffer.data(target + 1135);
    auto *t_1136 = buffer.data(target + 1136);
    auto *t_1137 = buffer.data(target + 1137);
    auto *t_1138 = buffer.data(target + 1138);
    auto *t_1139 = buffer.data(target + 1139);
    auto *t_1140 = buffer.data(target + 1140);
    auto *t_1141 = buffer.data(target + 1141);
    auto *t_1142 = buffer.data(target + 1142);
    auto *t_1143 = buffer.data(target + 1143);
    auto *t_1144 = buffer.data(target + 1144);
    auto *t_1145 = buffer.data(target + 1145);
    auto *t_1146 = buffer.data(target + 1146);
    auto *t_1147 = buffer.data(target + 1147);
    auto *t_1148 = buffer.data(target + 1148);
    auto *t_1149 = buffer.data(target + 1149);
    auto *t_1150 = buffer.data(target + 1150);
    auto *t_1151 = buffer.data(target + 1151);
    auto *t_1152 = buffer.data(target + 1152);
    auto *t_1153 = buffer.data(target + 1153);
    auto *t_1154 = buffer.data(target + 1154);
    auto *t_1155 = buffer.data(target + 1155);
    auto *t_1156 = buffer.data(target + 1156);
    auto *t_1157 = buffer.data(target + 1157);
    auto *t_1158 = buffer.data(target + 1158);
    auto *t_1159 = buffer.data(target + 1159);
    auto *t_1160 = buffer.data(target + 1160);
    auto *t_1161 = buffer.data(target + 1161);
    auto *t_1162 = buffer.data(target + 1162);
    auto *t_1163 = buffer.data(target + 1163);
    auto *t_1164 = buffer.data(target + 1164);
    auto *t_1165 = buffer.data(target + 1165);
    auto *t_1166 = buffer.data(target + 1166);
    auto *t_1167 = buffer.data(target + 1167);
    auto *t_1168 = buffer.data(target + 1168);
    auto *t_1169 = buffer.data(target + 1169);

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nsg0_975 = buffer.data(nsg0 + 975);
    const auto *nsg0_977 = buffer.data(nsg0 + 977);
    const auto *nsg0_980 = buffer.data(nsg0 + 980);
    const auto *nsg0_985 = buffer.data(nsg0 + 985);
    const auto *nsg0_987 = buffer.data(nsg0 + 987);
    const auto *nsg0_989 = buffer.data(nsg0 + 989);

    const auto *nsf_606 = buffer.data(nsf + 606);
    const auto *nsf_609 = buffer.data(nsf + 609);
    const auto *nsf_616 = buffer.data(nsf + 616);
    const auto *nsf_618 = buffer.data(nsf + 618);
    const auto *nsf_619 = buffer.data(nsf + 619);
    const auto *nsf_626 = buffer.data(nsf + 626);
    const auto *nsf_628 = buffer.data(nsf + 628);
    const auto *nsf_629 = buffer.data(nsf + 629);
    const auto *nsf_636 = buffer.data(nsf + 636);
    const auto *nsf_638 = buffer.data(nsf + 638);
    const auto *nsf_639 = buffer.data(nsf + 639);
    const auto *nsf_646 = buffer.data(nsf + 646);
    const auto *nsf_648 = buffer.data(nsf + 648);
    const auto *nsf_649 = buffer.data(nsf + 649);
    const auto *nsf_656 = buffer.data(nsf + 656);
    const auto *nsf_658 = buffer.data(nsf + 658);
    const auto *nsf_659 = buffer.data(nsf + 659);

    const auto *nsg1_975 = buffer.data(nsg1 + 975);
    const auto *nsg1_977 = buffer.data(nsg1 + 977);
    const auto *nsg1_980 = buffer.data(nsg1 + 980);
    const auto *nsg1_985 = buffer.data(nsg1 + 985);
    const auto *nsg1_987 = buffer.data(nsg1 + 987);
    const auto *nsg1_989 = buffer.data(nsg1 + 989);

    const auto *osd0_435 = buffer.data(osd0 + 435);
    const auto *osd0_437 = buffer.data(osd0 + 437);
    const auto *osd0_438 = buffer.data(osd0 + 438);
    const auto *osd0_439 = buffer.data(osd0 + 439);
    const auto *osd0_440 = buffer.data(osd0 + 440);
    const auto *osd0_441 = buffer.data(osd0 + 441);
    const auto *osd0_442 = buffer.data(osd0 + 442);
    const auto *osd0_443 = buffer.data(osd0 + 443);
    const auto *osd0_444 = buffer.data(osd0 + 444);
    const auto *osd0_445 = buffer.data(osd0 + 445);
    const auto *osd0_446 = buffer.data(osd0 + 446);
    const auto *osd0_447 = buffer.data(osd0 + 447);
    const auto *osd0_448 = buffer.data(osd0 + 448);
    const auto *osd0_449 = buffer.data(osd0 + 449);
    const auto *osd0_450 = buffer.data(osd0 + 450);
    const auto *osd0_451 = buffer.data(osd0 + 451);
    const auto *osd0_452 = buffer.data(osd0 + 452);
    const auto *osd0_453 = buffer.data(osd0 + 453);
    const auto *osd0_454 = buffer.data(osd0 + 454);
    const auto *osd0_455 = buffer.data(osd0 + 455);
    const auto *osd0_457 = buffer.data(osd0 + 457);
    const auto *osd0_459 = buffer.data(osd0 + 459);
    const auto *osd0_460 = buffer.data(osd0 + 460);
    const auto *osd0_462 = buffer.data(osd0 + 462);
    const auto *osd0_464 = buffer.data(osd0 + 464);
    const auto *osd0_465 = buffer.data(osd0 + 465);
    const auto *osd0_466 = buffer.data(osd0 + 466);
    const auto *osd0_467 = buffer.data(osd0 + 467);

    const auto *osd1_435 = buffer.data(osd1 + 435);
    const auto *osd1_437 = buffer.data(osd1 + 437);
    const auto *osd1_438 = buffer.data(osd1 + 438);
    const auto *osd1_439 = buffer.data(osd1 + 439);
    const auto *osd1_440 = buffer.data(osd1 + 440);
    const auto *osd1_441 = buffer.data(osd1 + 441);
    const auto *osd1_442 = buffer.data(osd1 + 442);
    const auto *osd1_443 = buffer.data(osd1 + 443);
    const auto *osd1_444 = buffer.data(osd1 + 444);
    const auto *osd1_445 = buffer.data(osd1 + 445);
    const auto *osd1_446 = buffer.data(osd1 + 446);
    const auto *osd1_447 = buffer.data(osd1 + 447);
    const auto *osd1_448 = buffer.data(osd1 + 448);
    const auto *osd1_449 = buffer.data(osd1 + 449);
    const auto *osd1_450 = buffer.data(osd1 + 450);
    const auto *osd1_451 = buffer.data(osd1 + 451);
    const auto *osd1_452 = buffer.data(osd1 + 452);
    const auto *osd1_453 = buffer.data(osd1 + 453);
    const auto *osd1_454 = buffer.data(osd1 + 454);
    const auto *osd1_455 = buffer.data(osd1 + 455);
    const auto *osd1_457 = buffer.data(osd1 + 457);
    const auto *osd1_459 = buffer.data(osd1 + 459);
    const auto *osd1_460 = buffer.data(osd1 + 460);
    const auto *osd1_462 = buffer.data(osd1 + 462);
    const auto *osd1_464 = buffer.data(osd1 + 464);
    const auto *osd1_465 = buffer.data(osd1 + 465);
    const auto *osd1_466 = buffer.data(osd1 + 466);
    const auto *osd1_467 = buffer.data(osd1 + 467);

    const auto *osf_726 = buffer.data(osf + 726);
    const auto *osf_727 = buffer.data(osf + 727);
    const auto *osf_728 = buffer.data(osf + 728);
    const auto *osf_729 = buffer.data(osf + 729);
    const auto *osf_730 = buffer.data(osf + 730);
    const auto *osf_731 = buffer.data(osf + 731);
    const auto *osf_732 = buffer.data(osf + 732);
    const auto *osf_733 = buffer.data(osf + 733);
    const auto *osf_734 = buffer.data(osf + 734);
    const auto *osf_735 = buffer.data(osf + 735);
    const auto *osf_736 = buffer.data(osf + 736);
    const auto *osf_737 = buffer.data(osf + 737);
    const auto *osf_738 = buffer.data(osf + 738);
    const auto *osf_739 = buffer.data(osf + 739);
    const auto *osf_740 = buffer.data(osf + 740);
    const auto *osf_741 = buffer.data(osf + 741);
    const auto *osf_742 = buffer.data(osf + 742);
    const auto *osf_743 = buffer.data(osf + 743);
    const auto *osf_744 = buffer.data(osf + 744);
    const auto *osf_745 = buffer.data(osf + 745);
    const auto *osf_746 = buffer.data(osf + 746);
    const auto *osf_747 = buffer.data(osf + 747);
    const auto *osf_748 = buffer.data(osf + 748);
    const auto *osf_749 = buffer.data(osf + 749);
    const auto *osf_750 = buffer.data(osf + 750);
    const auto *osf_751 = buffer.data(osf + 751);
    const auto *osf_752 = buffer.data(osf + 752);
    const auto *osf_753 = buffer.data(osf + 753);
    const auto *osf_754 = buffer.data(osf + 754);
    const auto *osf_755 = buffer.data(osf + 755);
    const auto *osf_756 = buffer.data(osf + 756);
    const auto *osf_757 = buffer.data(osf + 757);
    const auto *osf_758 = buffer.data(osf + 758);
    const auto *osf_759 = buffer.data(osf + 759);
    const auto *osf_761 = buffer.data(osf + 761);
    const auto *osf_763 = buffer.data(osf + 763);
    const auto *osf_764 = buffer.data(osf + 764);
    const auto *osf_766 = buffer.data(osf + 766);
    const auto *osf_767 = buffer.data(osf + 767);
    const auto *osf_768 = buffer.data(osf + 768);
    const auto *osf_769 = buffer.data(osf + 769);
    const auto *osf_770 = buffer.data(osf + 770);
    const auto *osf_772 = buffer.data(osf + 772);
    const auto *osf_773 = buffer.data(osf + 773);
    const auto *osf_775 = buffer.data(osf + 775);
    const auto *osf_776 = buffer.data(osf + 776);
    const auto *osf_777 = buffer.data(osf + 777);
    const auto *osf_778 = buffer.data(osf + 778);
    const auto *osf_779 = buffer.data(osf + 779);

#pragma omp simd aligned(t_1087, t_1088, t_1089, t_1090, t_1091, pc_x, pc_y, pc_z, nsf_606, \
                         nsf_616, osd0_435, osd1_435, osf_726, osf_727, osf_728, \
                         osf_729 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1087[k] = f_3 * pc_x[k] * osf_727[k];

        t_1088[k] = f_3 * pc_x[k] * osf_728[k];

        t_1089[k] = f_3 * pc_x[k] * osf_729[k];

        t_1090[k] = f_18 * nsf_616[k]
                    + f_1 * osd0_435[k]
                    - f_2 * osd1_435[k]
                    + f_3 * pc_y[k] * osf_726[k];

        t_1091[k] = f_17 * nsf_606[k]
                    + f_3 * pc_z[k] * osf_726[k];
    }

#pragma omp simd aligned(t_1092, t_1093, t_1094, pc_y, pc_z, nsf_609, nsf_618, nsf_619, \
                         osd0_437, osd1_437, osf_728, osf_729 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1092[k] = f_18 * nsf_618[k]
                    + f_4 * osd0_437[k]
                    - f_5 * osd1_437[k]
                    + f_3 * pc_y[k] * osf_728[k];

        t_1093[k] = f_18 * nsf_619[k]
                    + f_3 * pc_y[k] * osf_729[k];

        t_1094[k] = f_17 * nsf_609[k]
                    + f_1 * osd0_437[k]
                    - f_2 * osd1_437[k]
                    + f_3 * pc_z[k] * osf_729[k];
    }

#pragma omp simd aligned(t_1095, t_1096, t_1097, pc_x, osd0_438, osd0_439, osd0_440, osd1_438, \
                         osd1_439, osd1_440, osf_730, osf_731, \
                         osf_732 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1095[k] = f_1 * osd0_438[k]
                    - f_2 * osd1_438[k]
                    + f_3 * pc_x[k] * osf_730[k];

        t_1096[k] = f_10 * osd0_439[k]
                    - f_11 * osd1_439[k]
                    + f_3 * pc_x[k] * osf_731[k];

        t_1097[k] = f_10 * osd0_440[k]
                    - f_11 * osd1_440[k]
                    + f_3 * pc_x[k] * osf_732[k];
    }

#pragma omp simd aligned(t_1098, t_1099, t_1100, t_1101, pc_x, osd0_441, osd0_442, osd0_443, \
                         osd1_441, osd1_442, osd1_443, osf_733, osf_734, osf_735, \
                         osf_736 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1098[k] = f_4 * osd0_441[k]
                    - f_5 * osd1_441[k]
                    + f_3 * pc_x[k] * osf_733[k];

        t_1099[k] = f_4 * osd0_442[k]
                    - f_5 * osd1_442[k]
                    + f_3 * pc_x[k] * osf_734[k];

        t_1100[k] = f_4 * osd0_443[k]
                    - f_5 * osd1_443[k]
                    + f_3 * pc_x[k] * osf_735[k];

        t_1101[k] = f_3 * pc_x[k] * osf_736[k];
    }

#pragma omp simd aligned(t_1102, t_1103, t_1104, t_1105, t_1106, pc_x, pc_y, pc_z, nsf_616, \
                         nsf_626, osd0_441, osd1_441, osf_736, osf_737, osf_738, \
                         osf_739 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1102[k] = f_3 * pc_x[k] * osf_737[k];

        t_1103[k] = f_3 * pc_x[k] * osf_738[k];

        t_1104[k] = f_3 * pc_x[k] * osf_739[k];

        t_1105[k] = f_16 * nsf_626[k]
                    + f_1 * osd0_441[k]
                    - f_2 * osd1_441[k]
                    + f_3 * pc_y[k] * osf_736[k];

        t_1106[k] = f_15 * nsf_616[k]
                    + f_3 * pc_z[k] * osf_736[k];
    }

#pragma omp simd aligned(t_1107, t_1108, t_1109, pc_y, pc_z, nsf_619, nsf_628, nsf_629, \
                         osd0_443, osd1_443, osf_738, osf_739 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1107[k] = f_16 * nsf_628[k]
                    + f_4 * osd0_443[k]
                    - f_5 * osd1_443[k]
                    + f_3 * pc_y[k] * osf_738[k];

        t_1108[k] = f_16 * nsf_629[k]
                    + f_3 * pc_y[k] * osf_739[k];

        t_1109[k] = f_15 * nsf_619[k]
                    + f_1 * osd0_443[k]
                    - f_2 * osd1_443[k]
                    + f_3 * pc_z[k] * osf_739[k];
    }

#pragma omp simd aligned(t_1110, t_1111, t_1112, pc_x, osd0_444, osd0_445, osd0_446, osd1_444, \
                         osd1_445, osd1_446, osf_740, osf_741, \
                         osf_742 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1110[k] = f_1 * osd0_444[k]
                    - f_2 * osd1_444[k]
                    + f_3 * pc_x[k] * osf_740[k];

        t_1111[k] = f_10 * osd0_445[k]
                    - f_11 * osd1_445[k]
                    + f_3 * pc_x[k] * osf_741[k];

        t_1112[k] = f_10 * osd0_446[k]
                    - f_11 * osd1_446[k]
                    + f_3 * pc_x[k] * osf_742[k];
    }

#pragma omp simd aligned(t_1113, t_1114, t_1115, t_1116, pc_x, osd0_447, osd0_448, osd0_449, \
                         osd1_447, osd1_448, osd1_449, osf_743, osf_744, osf_745, \
                         osf_746 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1113[k] = f_4 * osd0_447[k]
                    - f_5 * osd1_447[k]
                    + f_3 * pc_x[k] * osf_743[k];

        t_1114[k] = f_4 * osd0_448[k]
                    - f_5 * osd1_448[k]
                    + f_3 * pc_x[k] * osf_744[k];

        t_1115[k] = f_4 * osd0_449[k]
                    - f_5 * osd1_449[k]
                    + f_3 * pc_x[k] * osf_745[k];

        t_1116[k] = f_3 * pc_x[k] * osf_746[k];
    }

#pragma omp simd aligned(t_1117, t_1118, t_1119, t_1120, t_1121, pc_x, pc_y, pc_z, nsf_626, \
                         nsf_636, osd0_447, osd1_447, osf_746, osf_747, osf_748, \
                         osf_749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1117[k] = f_3 * pc_x[k] * osf_747[k];

        t_1118[k] = f_3 * pc_x[k] * osf_748[k];

        t_1119[k] = f_3 * pc_x[k] * osf_749[k];

        t_1120[k] = f_14 * nsf_636[k]
                    + f_1 * osd0_447[k]
                    - f_2 * osd1_447[k]
                    + f_3 * pc_y[k] * osf_746[k];

        t_1121[k] = f_13 * nsf_626[k]
                    + f_3 * pc_z[k] * osf_746[k];
    }

#pragma omp simd aligned(t_1122, t_1123, t_1124, pc_y, pc_z, nsf_629, nsf_638, nsf_639, \
                         osd0_449, osd1_449, osf_748, osf_749 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1122[k] = f_14 * nsf_638[k]
                    + f_4 * osd0_449[k]
                    - f_5 * osd1_449[k]
                    + f_3 * pc_y[k] * osf_748[k];

        t_1123[k] = f_14 * nsf_639[k]
                    + f_3 * pc_y[k] * osf_749[k];

        t_1124[k] = f_13 * nsf_629[k]
                    + f_1 * osd0_449[k]
                    - f_2 * osd1_449[k]
                    + f_3 * pc_z[k] * osf_749[k];
    }

#pragma omp simd aligned(t_1125, t_1126, t_1127, pc_x, osd0_450, osd0_451, osd0_452, osd1_450, \
                         osd1_451, osd1_452, osf_750, osf_751, \
                         osf_752 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1125[k] = f_1 * osd0_450[k]
                    - f_2 * osd1_450[k]
                    + f_3 * pc_x[k] * osf_750[k];

        t_1126[k] = f_10 * osd0_451[k]
                    - f_11 * osd1_451[k]
                    + f_3 * pc_x[k] * osf_751[k];

        t_1127[k] = f_10 * osd0_452[k]
                    - f_11 * osd1_452[k]
                    + f_3 * pc_x[k] * osf_752[k];
    }

#pragma omp simd aligned(t_1128, t_1129, t_1130, t_1131, pc_x, osd0_453, osd0_454, osd0_455, \
                         osd1_453, osd1_454, osd1_455, osf_753, osf_754, osf_755, \
                         osf_756 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1128[k] = f_4 * osd0_453[k]
                    - f_5 * osd1_453[k]
                    + f_3 * pc_x[k] * osf_753[k];

        t_1129[k] = f_4 * osd0_454[k]
                    - f_5 * osd1_454[k]
                    + f_3 * pc_x[k] * osf_754[k];

        t_1130[k] = f_4 * osd0_455[k]
                    - f_5 * osd1_455[k]
                    + f_3 * pc_x[k] * osf_755[k];

        t_1131[k] = f_3 * pc_x[k] * osf_756[k];
    }

#pragma omp simd aligned(t_1132, t_1133, t_1134, t_1135, t_1136, pc_x, pc_y, pc_z, nsf_636, \
                         nsf_646, osd0_453, osd1_453, osf_756, osf_757, osf_758, \
                         osf_759 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1132[k] = f_3 * pc_x[k] * osf_757[k];

        t_1133[k] = f_3 * pc_x[k] * osf_758[k];

        t_1134[k] = f_3 * pc_x[k] * osf_759[k];

        t_1135[k] = f_8 * nsf_646[k]
                    + f_1 * osd0_453[k]
                    - f_2 * osd1_453[k]
                    + f_3 * pc_y[k] * osf_756[k];

        t_1136[k] = f_12 * nsf_636[k]
                    + f_3 * pc_z[k] * osf_756[k];
    }

#pragma omp simd aligned(t_1137, t_1138, t_1139, t_1140, pa_y, pc_y, pc_z, nsg0_975, nsf_639, \
                         nsf_648, nsf_649, nsg1_975, osd0_455, osd1_455, osf_758, \
                         osf_759 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1137[k] = f_8 * nsf_648[k]
                    + f_4 * osd0_455[k]
                    - f_5 * osd1_455[k]
                    + f_3 * pc_y[k] * osf_758[k];

        t_1138[k] = f_8 * nsf_649[k]
                    + f_3 * pc_y[k] * osf_759[k];

        t_1139[k] = f_12 * nsf_639[k]
                    + f_1 * osd0_455[k]
                    - f_2 * osd1_455[k]
                    + f_3 * pc_z[k] * osf_759[k];

        t_1140[k] = pa_y[k] * nsg0_975[k]
                    - f_6 * pc_y[k] * nsg1_975[k];
    }

#pragma omp simd aligned(t_1141, t_1142, t_1143, pa_y, pc_x, pc_y, nsg0_977, nsg1_977, \
                         osd0_457, osd0_459, osd1_457, osd1_459, osf_761, \
                         osf_763 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1141[k] = f_10 * osd0_457[k]
                    - f_11 * osd1_457[k]
                    + f_3 * pc_x[k] * osf_761[k];

        t_1142[k] = pa_y[k] * nsg0_977[k]
                    - f_6 * pc_y[k] * nsg1_977[k];

        t_1143[k] = f_4 * osd0_459[k]
                    - f_5 * osd1_459[k]
                    + f_3 * pc_x[k] * osf_763[k];
    }

#pragma omp simd aligned(t_1144, t_1145, t_1146, t_1147, t_1148, pa_y, pc_x, pc_y, nsg0_980, \
                         nsg1_980, osd0_460, osd1_460, osf_764, osf_766, osf_767, \
                         osf_768 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1144[k] = f_4 * osd0_460[k]
                    - f_5 * osd1_460[k]
                    + f_3 * pc_x[k] * osf_764[k];

        t_1145[k] = pa_y[k] * nsg0_980[k]
                    - f_6 * pc_y[k] * nsg1_980[k];

        t_1146[k] = f_3 * pc_x[k] * osf_766[k];

        t_1147[k] = f_3 * pc_x[k] * osf_767[k];

        t_1148[k] = f_3 * pc_x[k] * osf_768[k];
    }

#pragma omp simd aligned(t_1149, t_1150, t_1151, pa_y, pc_x, pc_y, pc_z, nsg0_985, nsf_646, \
                         nsf_656, nsg1_985, osf_766, osf_769 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1149[k] = f_3 * pc_x[k] * osf_769[k];

        t_1150[k] = pa_y[k] * nsg0_985[k]
                    + f_16 * nsf_656[k]
                    - f_6 * pc_y[k] * nsg1_985[k];

        t_1151[k] = f_9 * nsf_646[k]
                    + f_3 * pc_z[k] * osf_766[k];
    }

#pragma omp simd aligned(t_1152, t_1153, t_1154, pa_y, pc_y, nsg0_987, nsg0_989, nsf_658, \
                         nsf_659, nsg1_987, nsg1_989, osf_769 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1152[k] = pa_y[k] * nsg0_987[k]
                    + f_8 * nsf_658[k]
                    - f_6 * pc_y[k] * nsg1_987[k];

        t_1153[k] = f_7 * nsf_659[k]
                    + f_3 * pc_y[k] * osf_769[k];

        t_1154[k] = pa_y[k] * nsg0_989[k]
                    - f_6 * pc_y[k] * nsg1_989[k];
    }

#pragma omp simd aligned(t_1155, t_1156, t_1157, t_1158, t_1159, pc_x, pc_y, osd0_462, \
                         osd0_464, osd0_465, osd1_462, osd1_464, osd1_465, osf_770, osf_772, \
                         osf_773 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1155[k] = f_1 * osd0_462[k]
                    - f_2 * osd1_462[k]
                    + f_3 * pc_x[k] * osf_770[k];

        t_1156[k] = f_3 * pc_y[k] * osf_770[k];

        t_1157[k] = f_10 * osd0_464[k]
                    - f_11 * osd1_464[k]
                    + f_3 * pc_x[k] * osf_772[k];

        t_1158[k] = f_4 * osd0_465[k]
                    - f_5 * osd1_465[k]
                    + f_3 * pc_x[k] * osf_773[k];

        t_1159[k] = f_3 * pc_y[k] * osf_772[k];
    }

#pragma omp simd aligned(t_1160, t_1161, t_1162, t_1163, t_1164, pc_x, osd0_467, osd1_467, \
                         osf_775, osf_776, osf_777, osf_778, osf_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1160[k] = f_4 * osd0_467[k]
                    - f_5 * osd1_467[k]
                    + f_3 * pc_x[k] * osf_775[k];

        t_1161[k] = f_3 * pc_x[k] * osf_776[k];

        t_1162[k] = f_3 * pc_x[k] * osf_777[k];

        t_1163[k] = f_3 * pc_x[k] * osf_778[k];

        t_1164[k] = f_3 * pc_x[k] * osf_779[k];
    }

#pragma omp simd aligned(t_1165, t_1166, t_1167, t_1168, pc_y, osd0_465, osd0_466, osd0_467, \
                         osd1_465, osd1_466, osd1_467, osf_776, osf_777, osf_778, \
                         osf_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1165[k] = f_1 * osd0_465[k]
                    - f_2 * osd1_465[k]
                    + f_3 * pc_y[k] * osf_776[k];

        t_1166[k] = f_10 * osd0_466[k]
                    - f_11 * osd1_466[k]
                    + f_3 * pc_y[k] * osf_777[k];

        t_1167[k] = f_4 * osd0_467[k]
                    - f_5 * osd1_467[k]
                    + f_3 * pc_y[k] * osf_778[k];

        t_1168[k] = f_3 * pc_y[k] * osf_779[k];
    }

#pragma omp simd aligned(t_1169, pc_z, nsf_659, osd0_467, osd1_467, \
                         osf_779 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_1169[k] = f_0 * nsf_659[k]
                    + f_1 * osd0_467[k]
                    - f_2 * osd1_467[k]
                    + f_3 * pc_z[k] * osf_779[k];
    }
}

auto
compute_prim_osg_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t nsg0, const size_t nsf,
                                                   const size_t nsg1, const size_t osd0,
                                                   const size_t osd1, const size_t osf,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_osg_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, nsg0, nsf,
                                                              nsg1, osd0, osd1, osf, ncols,
                                                              gamma, p, q);

    compute_prim_osg_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, nsg0, nsf,
                                                              nsg1, osd0, osd1, osf, ncols,
                                                              gamma, p, q);

    compute_prim_osg_three_center_electron_repulsion_0_piece2(buffer, target, pa, pc, nsg0, nsf,
                                                              nsg1, osd0, osd1, osf, ncols,
                                                              gamma, p, q);

    compute_prim_osg_three_center_electron_repulsion_0_piece3(buffer, target, pa, pc, nsg0, nsf,
                                                              nsg1, osd0, osd1, osf, ncols,
                                                              gamma, p, q);

    compute_prim_osg_three_center_electron_repulsion_0_piece4(buffer, target, pa, pc, nsg0, nsf,
                                                              nsg1, osd0, osd1, osf, ncols,
                                                              gamma, p, q);

    compute_prim_osg_three_center_electron_repulsion_0_piece5(buffer, target, pa, pc, nsg0, nsf,
                                                              nsg1, osd0, osd1, osf, ncols,
                                                              gamma, p, q);

    compute_prim_osg_three_center_electron_repulsion_0_piece6(buffer, target, pa, pc, nsg0, nsf,
                                                              nsg1, osd0, osd1, osf, ncols,
                                                              gamma, p, q);

    compute_prim_osg_three_center_electron_repulsion_0_piece7(buffer, target, pa, pc, nsg0, nsf,
                                                              nsg1, osf, ncols, gamma, p, q);

    compute_prim_osg_three_center_electron_repulsion_0_piece8(buffer, target, pa, pc, nsg0, nsf,
                                                              nsg1, osd0, osd1, osf, ncols,
                                                              gamma, p, q);

    compute_prim_osg_three_center_electron_repulsion_0_piece9(buffer, target, pa, pc, nsg0, nsf,
                                                              nsg1, osd0, osd1, osf, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
