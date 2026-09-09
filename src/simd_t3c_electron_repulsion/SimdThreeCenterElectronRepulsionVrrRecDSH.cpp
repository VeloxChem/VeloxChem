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


#include "SimdThreeCenterElectronRepulsionVrrRecDSH.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_dsh_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t psh0, const size_t psg,
                                                   const size_t psh1, const size_t dsf0,
                                                   const size_t dsf1, const size_t dsg,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.5 / q;
    const auto f_11 = 1.5 / gamma;
    const auto f_12 = 1.5 * p / (gamma * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *psh0_0 = buffer.data(psh0 + 0);
    const auto *psh0_3 = buffer.data(psh0 + 3);
    const auto *psh0_5 = buffer.data(psh0 + 5);
    const auto *psh0_6 = buffer.data(psh0 + 6);
    const auto *psh0_9 = buffer.data(psh0 + 9);
    const auto *psh0_22 = buffer.data(psh0 + 22);
    const auto *psh0_24 = buffer.data(psh0 + 24);
    const auto *psh0_27 = buffer.data(psh0 + 27);
    const auto *psh0_36 = buffer.data(psh0 + 36);
    const auto *psh0_38 = buffer.data(psh0 + 38);
    const auto *psh0_39 = buffer.data(psh0 + 39);
    const auto *psh0_41 = buffer.data(psh0 + 41);
    const auto *psh0_42 = buffer.data(psh0 + 42);
    const auto *psh0_44 = buffer.data(psh0 + 44);
    const auto *psh0_47 = buffer.data(psh0 + 47);
    const auto *psh0_51 = buffer.data(psh0 + 51);
    const auto *psh0_57 = buffer.data(psh0 + 57);
    const auto *psh0_58 = buffer.data(psh0 + 58);
    const auto *psh0_59 = buffer.data(psh0 + 59);
    const auto *psh0_60 = buffer.data(psh0 + 60);
    const auto *psh0_62 = buffer.data(psh0 + 62);

    const auto *psg_0 = buffer.data(psg + 0);
    const auto *psg_5 = buffer.data(psg + 5);
    const auto *psg_10 = buffer.data(psg + 10);
    const auto *psg_12 = buffer.data(psg + 12);
    const auto *psg_14 = buffer.data(psg + 14);
    const auto *psg_18 = buffer.data(psg + 18);
    const auto *psg_21 = buffer.data(psg + 21);
    const auto *psg_25 = buffer.data(psg + 25);
    const auto *psg_27 = buffer.data(psg + 27);
    const auto *psg_28 = buffer.data(psg + 28);
    const auto *psg_29 = buffer.data(psg + 29);
    const auto *psg_35 = buffer.data(psg + 35);
    const auto *psg_39 = buffer.data(psg + 39);
    const auto *psg_40 = buffer.data(psg + 40);
    const auto *psg_41 = buffer.data(psg + 41);
    const auto *psg_42 = buffer.data(psg + 42);
    const auto *psg_43 = buffer.data(psg + 43);
    const auto *psg_44 = buffer.data(psg + 44);

    const auto *psh1_0 = buffer.data(psh1 + 0);
    const auto *psh1_3 = buffer.data(psh1 + 3);
    const auto *psh1_5 = buffer.data(psh1 + 5);
    const auto *psh1_6 = buffer.data(psh1 + 6);
    const auto *psh1_9 = buffer.data(psh1 + 9);
    const auto *psh1_22 = buffer.data(psh1 + 22);
    const auto *psh1_24 = buffer.data(psh1 + 24);
    const auto *psh1_27 = buffer.data(psh1 + 27);
    const auto *psh1_36 = buffer.data(psh1 + 36);
    const auto *psh1_38 = buffer.data(psh1 + 38);
    const auto *psh1_39 = buffer.data(psh1 + 39);
    const auto *psh1_41 = buffer.data(psh1 + 41);
    const auto *psh1_42 = buffer.data(psh1 + 42);
    const auto *psh1_44 = buffer.data(psh1 + 44);
    const auto *psh1_47 = buffer.data(psh1 + 47);
    const auto *psh1_51 = buffer.data(psh1 + 51);
    const auto *psh1_57 = buffer.data(psh1 + 57);
    const auto *psh1_58 = buffer.data(psh1 + 58);
    const auto *psh1_59 = buffer.data(psh1 + 59);
    const auto *psh1_60 = buffer.data(psh1 + 60);
    const auto *psh1_62 = buffer.data(psh1 + 62);

    const auto *dsf0_0 = buffer.data(dsf0 + 0);
    const auto *dsf0_1 = buffer.data(dsf0 + 1);
    const auto *dsf0_2 = buffer.data(dsf0 + 2);
    const auto *dsf0_6 = buffer.data(dsf0 + 6);
    const auto *dsf0_8 = buffer.data(dsf0 + 8);
    const auto *dsf0_9 = buffer.data(dsf0 + 9);
    const auto *dsf0_22 = buffer.data(dsf0 + 22);
    const auto *dsf0_30 = buffer.data(dsf0 + 30);
    const auto *dsf0_31 = buffer.data(dsf0 + 31);
    const auto *dsf0_33 = buffer.data(dsf0 + 33);
    const auto *dsf0_35 = buffer.data(dsf0 + 35);
    const auto *dsf0_36 = buffer.data(dsf0 + 36);
    const auto *dsf0_37 = buffer.data(dsf0 + 37);
    const auto *dsf0_38 = buffer.data(dsf0 + 38);
    const auto *dsf0_39 = buffer.data(dsf0 + 39);
    const auto *dsf0_44 = buffer.data(dsf0 + 44);
    const auto *dsf0_47 = buffer.data(dsf0 + 47);
    const auto *dsf0_48 = buffer.data(dsf0 + 48);
    const auto *dsf0_50 = buffer.data(dsf0 + 50);
    const auto *dsf0_52 = buffer.data(dsf0 + 52);
    const auto *dsf0_53 = buffer.data(dsf0 + 53);
    const auto *dsf0_55 = buffer.data(dsf0 + 55);
    const auto *dsf0_56 = buffer.data(dsf0 + 56);
    const auto *dsf0_57 = buffer.data(dsf0 + 57);
    const auto *dsf0_58 = buffer.data(dsf0 + 58);
    const auto *dsf0_59 = buffer.data(dsf0 + 59);

    const auto *dsf1_0 = buffer.data(dsf1 + 0);
    const auto *dsf1_1 = buffer.data(dsf1 + 1);
    const auto *dsf1_2 = buffer.data(dsf1 + 2);
    const auto *dsf1_6 = buffer.data(dsf1 + 6);
    const auto *dsf1_8 = buffer.data(dsf1 + 8);
    const auto *dsf1_9 = buffer.data(dsf1 + 9);
    const auto *dsf1_22 = buffer.data(dsf1 + 22);
    const auto *dsf1_30 = buffer.data(dsf1 + 30);
    const auto *dsf1_31 = buffer.data(dsf1 + 31);
    const auto *dsf1_33 = buffer.data(dsf1 + 33);
    const auto *dsf1_35 = buffer.data(dsf1 + 35);
    const auto *dsf1_36 = buffer.data(dsf1 + 36);
    const auto *dsf1_37 = buffer.data(dsf1 + 37);
    const auto *dsf1_38 = buffer.data(dsf1 + 38);
    const auto *dsf1_39 = buffer.data(dsf1 + 39);
    const auto *dsf1_44 = buffer.data(dsf1 + 44);
    const auto *dsf1_47 = buffer.data(dsf1 + 47);
    const auto *dsf1_48 = buffer.data(dsf1 + 48);
    const auto *dsf1_50 = buffer.data(dsf1 + 50);
    const auto *dsf1_52 = buffer.data(dsf1 + 52);
    const auto *dsf1_53 = buffer.data(dsf1 + 53);
    const auto *dsf1_55 = buffer.data(dsf1 + 55);
    const auto *dsf1_56 = buffer.data(dsf1 + 56);
    const auto *dsf1_57 = buffer.data(dsf1 + 57);
    const auto *dsf1_58 = buffer.data(dsf1 + 58);
    const auto *dsf1_59 = buffer.data(dsf1 + 59);

    const auto *dsg_0 = buffer.data(dsg + 0);
    const auto *dsg_1 = buffer.data(dsg + 1);
    const auto *dsg_2 = buffer.data(dsg + 2);
    const auto *dsg_3 = buffer.data(dsg + 3);
    const auto *dsg_5 = buffer.data(dsg + 5);
    const auto *dsg_6 = buffer.data(dsg + 6);
    const auto *dsg_9 = buffer.data(dsg + 9);
    const auto *dsg_10 = buffer.data(dsg + 10);
    const auto *dsg_12 = buffer.data(dsg + 12);
    const auto *dsg_13 = buffer.data(dsg + 13);
    const auto *dsg_14 = buffer.data(dsg + 14);
    const auto *dsg_15 = buffer.data(dsg + 15);
    const auto *dsg_16 = buffer.data(dsg + 16);
    const auto *dsg_18 = buffer.data(dsg + 18);
    const auto *dsg_20 = buffer.data(dsg + 20);
    const auto *dsg_21 = buffer.data(dsg + 21);
    const auto *dsg_25 = buffer.data(dsg + 25);
    const auto *dsg_27 = buffer.data(dsg + 27);
    const auto *dsg_28 = buffer.data(dsg + 28);
    const auto *dsg_29 = buffer.data(dsg + 29);
    const auto *dsg_30 = buffer.data(dsg + 30);
    const auto *dsg_32 = buffer.data(dsg + 32);
    const auto *dsg_34 = buffer.data(dsg + 34);
    const auto *dsg_35 = buffer.data(dsg + 35);
    const auto *dsg_39 = buffer.data(dsg + 39);
    const auto *dsg_40 = buffer.data(dsg + 40);
    const auto *dsg_41 = buffer.data(dsg + 41);
    const auto *dsg_42 = buffer.data(dsg + 42);
    const auto *dsg_44 = buffer.data(dsg + 44);
    const auto *dsg_45 = buffer.data(dsg + 45);
    const auto *dsg_46 = buffer.data(dsg + 46);
    const auto *dsg_48 = buffer.data(dsg + 48);
    const auto *dsg_50 = buffer.data(dsg + 50);
    const auto *dsg_51 = buffer.data(dsg + 51);
    const auto *dsg_53 = buffer.data(dsg + 53);
    const auto *dsg_54 = buffer.data(dsg + 54);
    const auto *dsg_55 = buffer.data(dsg + 55);
    const auto *dsg_56 = buffer.data(dsg + 56);
    const auto *dsg_57 = buffer.data(dsg + 57);
    const auto *dsg_58 = buffer.data(dsg + 58);
    const auto *dsg_59 = buffer.data(dsg + 59);
    const auto *dsg_64 = buffer.data(dsg + 64);
    const auto *dsg_67 = buffer.data(dsg + 67);
    const auto *dsg_68 = buffer.data(dsg + 68);
    const auto *dsg_70 = buffer.data(dsg + 70);
    const auto *dsg_71 = buffer.data(dsg + 71);
    const auto *dsg_72 = buffer.data(dsg + 72);
    const auto *dsg_73 = buffer.data(dsg + 73);
    const auto *dsg_74 = buffer.data(dsg + 74);
    const auto *dsg_75 = buffer.data(dsg + 75);
    const auto *dsg_77 = buffer.data(dsg + 77);
    const auto *dsg_78 = buffer.data(dsg + 78);
    const auto *dsg_80 = buffer.data(dsg + 80);
    const auto *dsg_81 = buffer.data(dsg + 81);
    const auto *dsg_82 = buffer.data(dsg + 82);
    const auto *dsg_84 = buffer.data(dsg + 84);
    const auto *dsg_85 = buffer.data(dsg + 85);
    const auto *dsg_86 = buffer.data(dsg + 86);
    const auto *dsg_87 = buffer.data(dsg + 87);
    const auto *dsg_88 = buffer.data(dsg + 88);
    const auto *dsg_89 = buffer.data(dsg + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, psg_0, dsf0_0, \
                         dsf1_0, dsg_0, dsg_1, dsg_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * psg_0[k]
                 + f_1 * dsf0_0[k]
                 - f_2 * dsf1_0[k]
                 + f_3 * pc_x[k] * dsg_0[k];

        t_1[k] = f_3 * pc_y[k] * dsg_0[k];

        t_2[k] = f_3 * pc_z[k] * dsg_0[k];

        t_3[k] = f_4 * dsf0_0[k]
                 - f_5 * dsf1_0[k]
                 + f_3 * pc_y[k] * dsg_1[k];

        t_4[k] = f_3 * pc_y[k] * dsg_2[k];

        t_5[k] = f_4 * dsf0_0[k]
                 - f_5 * dsf1_0[k]
                 + f_3 * pc_z[k] * dsg_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_x, pc_y, pc_z, psg_10, dsf0_1, dsf0_2, \
                         dsf1_1, dsf1_2, dsg_3, dsg_5, dsg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_6 * dsf0_1[k]
                 - f_7 * dsf1_1[k]
                 + f_3 * pc_y[k] * dsg_3[k];

        t_7[k] = f_3 * pc_z[k] * dsg_3[k];

        t_8[k] = f_3 * pc_y[k] * dsg_5[k];

        t_9[k] = f_6 * dsf0_2[k]
                 - f_7 * dsf1_2[k]
                 + f_3 * pc_z[k] * dsg_5[k];

        t_10[k] = f_0 * psg_10[k]
                  + f_3 * pc_x[k] * dsg_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pc_x, pc_y, pc_z, psg_12, psg_14, dsg_6, \
                         dsg_9, dsg_12, dsg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_3 * pc_z[k] * dsg_6[k];

        t_12[k] = f_0 * psg_12[k]
                  + f_3 * pc_x[k] * dsg_12[k];

        t_13[k] = f_3 * pc_y[k] * dsg_9[k];

        t_14[k] = f_0 * psg_14[k]
                  + f_3 * pc_x[k] * dsg_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pc_y, pc_z, dsf0_6, dsf0_8, dsf0_9, dsf1_6, \
                         dsf1_8, dsf1_9, dsg_10, dsg_12, dsg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * dsf0_6[k]
                  - f_2 * dsf1_6[k]
                  + f_3 * pc_y[k] * dsg_10[k];

        t_16[k] = f_3 * pc_z[k] * dsg_10[k];

        t_17[k] = f_6 * dsf0_8[k]
                  - f_7 * dsf1_8[k]
                  + f_3 * pc_y[k] * dsg_12[k];

        t_18[k] = f_4 * dsf0_9[k]
                  - f_5 * dsf1_9[k]
                  + f_3 * pc_y[k] * dsg_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pa_y, pc_y, pc_z, psh0_0, psg_0, \
                         psh1_0, dsf0_9, dsf1_9, dsg_14, dsg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * dsg_14[k];

        t_20[k] = f_1 * dsf0_9[k]
                  - f_2 * dsf1_9[k]
                  + f_3 * pc_z[k] * dsg_14[k];

        t_21[k] = pa_y[k] * psh0_0[k]
                  - f_8 * pc_y[k] * psh1_0[k];

        t_22[k] = f_9 * psg_0[k]
                  + f_3 * pc_y[k] * dsg_15[k];

        t_23[k] = f_3 * pc_z[k] * dsg_15[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pa_x, pa_y, pc_x, pc_y, pc_z, psh0_5, psh0_24, \
                         psg_18, psh1_5, psh1_24, dsg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pa_x[k] * psh0_24[k]
                  + f_10 * psg_18[k]
                  - f_8 * pc_x[k] * psh1_24[k];

        t_25[k] = f_3 * pc_z[k] * dsg_16[k];

        t_26[k] = pa_y[k] * psh0_5[k]
                  - f_8 * pc_y[k] * psh1_5[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_x, pc_x, pc_y, pc_z, psh0_27, psg_5, psg_21, \
                         psh1_27, dsg_18, dsg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pa_x[k] * psh0_27[k]
                  + f_0 * psg_21[k]
                  - f_8 * pc_x[k] * psh1_27[k];

        t_28[k] = f_3 * pc_z[k] * dsg_18[k];

        t_29[k] = f_9 * psg_5[k]
                  + f_3 * pc_y[k] * dsg_20[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_y, pc_x, pc_y, pc_z, psh0_9, psg_25, \
                         psg_27, psh1_9, dsg_21, dsg_25, dsg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pa_y[k] * psh0_9[k]
                  - f_8 * pc_y[k] * psh1_9[k];

        t_31[k] = f_9 * psg_25[k]
                  + f_3 * pc_x[k] * dsg_25[k];

        t_32[k] = f_3 * pc_z[k] * dsg_21[k];

        t_33[k] = f_9 * psg_27[k]
                  + f_3 * pc_x[k] * dsg_27[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pa_x, pc_x, pc_z, psh0_36, psg_28, psg_29, \
                         psh1_36, dsg_25, dsg_28, dsg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_9 * psg_28[k]
                  + f_3 * pc_x[k] * dsg_28[k];

        t_35[k] = f_9 * psg_29[k]
                  + f_3 * pc_x[k] * dsg_29[k];

        t_36[k] = pa_x[k] * psh0_36[k]
                  - f_8 * pc_x[k] * psh1_36[k];

        t_37[k] = f_3 * pc_z[k] * dsg_25[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_x, pc_x, pc_y, psh0_38, psh0_39, psh0_41, \
                         psg_14, psh1_38, psh1_39, psh1_41, dsg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pa_x[k] * psh0_38[k]
                  - f_8 * pc_x[k] * psh1_38[k];

        t_39[k] = pa_x[k] * psh0_39[k]
                  - f_8 * pc_x[k] * psh1_39[k];

        t_40[k] = f_9 * psg_14[k]
                  + f_3 * pc_y[k] * dsg_29[k];

        t_41[k] = pa_x[k] * psh0_41[k]
                  - f_8 * pc_x[k] * psh1_41[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pa_z, pc_y, pc_z, psh0_0, psh0_3, \
                         psg_0, psh1_0, psh1_3, dsg_30, dsg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_z[k] * psh0_0[k]
                  - f_8 * pc_z[k] * psh1_0[k];

        t_43[k] = f_3 * pc_y[k] * dsg_30[k];

        t_44[k] = f_9 * psg_0[k]
                  + f_3 * pc_z[k] * dsg_30[k];

        t_45[k] = pa_z[k] * psh0_3[k]
                  - f_8 * pc_z[k] * psh1_3[k];

        t_46[k] = f_3 * pc_y[k] * dsg_32[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pa_x, pa_z, pc_x, pc_y, pc_z, psh0_6, psh0_47, \
                         psg_35, psh1_6, psh1_47, dsf0_22, dsf1_22, \
                         dsg_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = pa_x[k] * psh0_47[k]
                  + f_10 * psg_35[k]
                  - f_8 * pc_x[k] * psh1_47[k];

        t_48[k] = pa_z[k] * psh0_6[k]
                  - f_8 * pc_z[k] * psh1_6[k];

        t_49[k] = f_4 * dsf0_22[k]
                  - f_5 * dsf1_22[k]
                  + f_3 * pc_y[k] * dsg_34[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_x, pc_x, pc_y, psh0_51, psg_39, psg_40, \
                         psg_41, psh1_51, dsg_35, dsg_40, dsg_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_3 * pc_y[k] * dsg_35[k];

        t_51[k] = pa_x[k] * psh0_51[k]
                  + f_0 * psg_39[k]
                  - f_8 * pc_x[k] * psh1_51[k];

        t_52[k] = f_9 * psg_40[k]
                  + f_3 * pc_x[k] * dsg_40[k];

        t_53[k] = f_9 * psg_41[k]
                  + f_3 * pc_x[k] * dsg_41[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_x, pc_x, pc_y, psh0_57, psg_42, psg_44, \
                         psh1_57, dsg_39, dsg_42, dsg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_9 * psg_42[k]
                  + f_3 * pc_x[k] * dsg_42[k];

        t_55[k] = f_3 * pc_y[k] * dsg_39[k];

        t_56[k] = f_9 * psg_44[k]
                  + f_3 * pc_x[k] * dsg_44[k];

        t_57[k] = pa_x[k] * psh0_57[k]
                  - f_8 * pc_x[k] * psh1_57[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pa_x, pc_x, pc_y, psh0_58, psh0_59, psh0_60, \
                         psh1_58, psh1_59, psh1_60, dsg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = pa_x[k] * psh0_58[k]
                  - f_8 * pc_x[k] * psh1_58[k];

        t_59[k] = pa_x[k] * psh0_59[k]
                  - f_8 * pc_x[k] * psh1_59[k];

        t_60[k] = pa_x[k] * psh0_60[k]
                  - f_8 * pc_x[k] * psh1_60[k];

        t_61[k] = f_3 * pc_y[k] * dsg_44[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_x, pc_x, pc_z, psh0_62, psh1_62, dsf0_30, \
                         dsf0_31, dsf1_30, dsf1_31, dsg_45, dsg_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pa_x[k] * psh0_62[k]
                  - f_8 * pc_x[k] * psh1_62[k];

        t_63[k] = f_1 * dsf0_30[k]
                  - f_2 * dsf1_30[k]
                  + f_3 * pc_x[k] * dsg_45[k];

        t_64[k] = f_11 * dsf0_31[k]
                  - f_12 * dsf1_31[k]
                  + f_3 * pc_x[k] * dsg_46[k];

        t_65[k] = f_3 * pc_z[k] * dsg_45[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pc_x, pc_z, dsf0_33, dsf0_35, dsf0_36, \
                         dsf1_33, dsf1_35, dsf1_36, dsg_46, dsg_48, dsg_50, \
                         dsg_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_6 * dsf0_33[k]
                  - f_7 * dsf1_33[k]
                  + f_3 * pc_x[k] * dsg_48[k];

        t_67[k] = f_3 * pc_z[k] * dsg_46[k];

        t_68[k] = f_6 * dsf0_35[k]
                  - f_7 * dsf1_35[k]
                  + f_3 * pc_x[k] * dsg_50[k];

        t_69[k] = f_4 * dsf0_36[k]
                  - f_5 * dsf1_36[k]
                  + f_3 * pc_x[k] * dsg_51[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, t_74, pc_x, pc_z, dsf0_38, dsf0_39, dsf1_38, \
                         dsf1_39, dsg_48, dsg_53, dsg_54, dsg_55, \
                         dsg_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_3 * pc_z[k] * dsg_48[k];

        t_71[k] = f_4 * dsf0_38[k]
                  - f_5 * dsf1_38[k]
                  + f_3 * pc_x[k] * dsg_53[k];

        t_72[k] = f_4 * dsf0_39[k]
                  - f_5 * dsf1_39[k]
                  + f_3 * pc_x[k] * dsg_54[k];

        t_73[k] = f_3 * pc_x[k] * dsg_55[k];

        t_74[k] = f_3 * pc_x[k] * dsg_56[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, pc_x, pc_y, pc_z, psg_25, dsf0_36, \
                         dsf1_36, dsg_55, dsg_57, dsg_58, dsg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_3 * pc_x[k] * dsg_57[k];

        t_76[k] = f_3 * pc_x[k] * dsg_58[k];

        t_77[k] = f_3 * pc_x[k] * dsg_59[k];

        t_78[k] = f_0 * psg_25[k]
                  + f_1 * dsf0_36[k]
                  - f_2 * dsf1_36[k]
                  + f_3 * pc_y[k] * dsg_55[k];

        t_79[k] = f_3 * pc_z[k] * dsg_55[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pc_y, pc_z, psg_29, dsf0_36, dsf0_37, \
                         dsf0_39, dsf1_36, dsf1_37, dsf1_39, dsg_56, dsg_57, \
                         dsg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_4 * dsf0_36[k]
                  - f_5 * dsf1_36[k]
                  + f_3 * pc_z[k] * dsg_56[k];

        t_81[k] = f_6 * dsf0_37[k]
                  - f_7 * dsf1_37[k]
                  + f_3 * pc_z[k] * dsg_57[k];

        t_82[k] = f_0 * psg_29[k]
                  + f_3 * pc_y[k] * dsg_59[k];

        t_83[k] = f_1 * dsf0_39[k]
                  - f_2 * dsf1_39[k]
                  + f_3 * pc_z[k] * dsg_59[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pa_y, pa_z, pc_y, pc_z, psh0_22, psh0_24, \
                         psh0_42, psh0_44, psh1_22, psh1_24, psh1_42, \
                         psh1_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = pa_y[k] * psh0_42[k]
                  - f_8 * pc_y[k] * psh1_42[k];

        t_85[k] = pa_z[k] * psh0_22[k]
                  - f_8 * pc_z[k] * psh1_22[k];

        t_86[k] = pa_y[k] * psh0_44[k]
                  - f_8 * pc_y[k] * psh1_44[k];

        t_87[k] = pa_z[k] * psh0_24[k]
                  - f_8 * pc_z[k] * psh1_24[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, pa_y, pa_z, pc_x, pc_y, pc_z, psh0_27, psh0_47, \
                         psh1_27, psh1_47, dsf0_44, dsf1_44, dsg_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_6 * dsf0_44[k]
                  - f_7 * dsf1_44[k]
                  + f_3 * pc_x[k] * dsg_64[k];

        t_89[k] = pa_y[k] * psh0_47[k]
                  - f_8 * pc_y[k] * psh1_47[k];

        t_90[k] = pa_z[k] * psh0_27[k]
                  - f_8 * pc_z[k] * psh1_27[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pa_y, pc_x, pc_y, psh0_51, psh1_51, dsf0_47, \
                         dsf0_48, dsf1_47, dsf1_48, dsg_67, dsg_68, \
                         dsg_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_4 * dsf0_47[k]
                  - f_5 * dsf1_47[k]
                  + f_3 * pc_x[k] * dsg_67[k];

        t_92[k] = f_4 * dsf0_48[k]
                  - f_5 * dsf1_48[k]
                  + f_3 * pc_x[k] * dsg_68[k];

        t_93[k] = pa_y[k] * psh0_51[k]
                  - f_8 * pc_y[k] * psh1_51[k];

        t_94[k] = f_3 * pc_x[k] * dsg_70[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, pa_z, pc_x, pc_z, psh0_36, psh1_36, \
                         dsg_71, dsg_72, dsg_73, dsg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_3 * pc_x[k] * dsg_71[k];

        t_96[k] = f_3 * pc_x[k] * dsg_72[k];

        t_97[k] = f_3 * pc_x[k] * dsg_73[k];

        t_98[k] = f_3 * pc_x[k] * dsg_74[k];

        t_99[k] = pa_z[k] * psh0_36[k]
                  - f_8 * pc_z[k] * psh1_36[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, pa_y, pc_y, pc_z, psh0_59, psh0_60, psg_25, \
                         psg_42, psg_43, psh1_59, psh1_60, dsg_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_9 * psg_25[k]
                   + f_3 * pc_z[k] * dsg_70[k];

        t_101[k] = pa_y[k] * psh0_59[k]
                   + f_10 * psg_42[k]
                   - f_8 * pc_y[k] * psh1_59[k];

        t_102[k] = pa_y[k] * psh0_60[k]
                   + f_0 * psg_43[k]
                   - f_8 * pc_y[k] * psh1_60[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, pa_y, pc_x, pc_y, psh0_62, psg_44, \
                         psh1_62, dsf0_50, dsf1_50, dsg_74, dsg_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_9 * psg_44[k]
                   + f_3 * pc_y[k] * dsg_74[k];

        t_104[k] = pa_y[k] * psh0_62[k]
                   - f_8 * pc_y[k] * psh1_62[k];

        t_105[k] = f_1 * dsf0_50[k]
                   - f_2 * dsf1_50[k]
                   + f_3 * pc_x[k] * dsg_75[k];

        t_106[k] = f_3 * pc_y[k] * dsg_75[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pc_x, pc_y, dsf0_52, dsf0_53, dsf0_55, \
                         dsf1_52, dsf1_53, dsf1_55, dsg_77, dsg_78, \
                         dsg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_11 * dsf0_52[k]
                   - f_12 * dsf1_52[k]
                   + f_3 * pc_x[k] * dsg_77[k];

        t_108[k] = f_6 * dsf0_53[k]
                   - f_7 * dsf1_53[k]
                   + f_3 * pc_x[k] * dsg_78[k];

        t_109[k] = f_3 * pc_y[k] * dsg_77[k];

        t_110[k] = f_6 * dsf0_55[k]
                   - f_7 * dsf1_55[k]
                   + f_3 * pc_x[k] * dsg_80[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pc_x, pc_y, dsf0_56, dsf0_57, dsf0_59, \
                         dsf1_56, dsf1_57, dsf1_59, dsg_80, dsg_81, dsg_82, \
                         dsg_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_4 * dsf0_56[k]
                   - f_5 * dsf1_56[k]
                   + f_3 * pc_x[k] * dsg_81[k];

        t_112[k] = f_4 * dsf0_57[k]
                   - f_5 * dsf1_57[k]
                   + f_3 * pc_x[k] * dsg_82[k];

        t_113[k] = f_3 * pc_y[k] * dsg_80[k];

        t_114[k] = f_4 * dsf0_59[k]
                   - f_5 * dsf1_59[k]
                   + f_3 * pc_x[k] * dsg_84[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, t_120, pc_x, pc_y, dsf0_56, \
                         dsf1_56, dsg_85, dsg_86, dsg_87, dsg_88, \
                         dsg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_3 * pc_x[k] * dsg_85[k];

        t_116[k] = f_3 * pc_x[k] * dsg_86[k];

        t_117[k] = f_3 * pc_x[k] * dsg_87[k];

        t_118[k] = f_3 * pc_x[k] * dsg_88[k];

        t_119[k] = f_3 * pc_x[k] * dsg_89[k];

        t_120[k] = f_1 * dsf0_56[k]
                   - f_2 * dsf1_56[k]
                   + f_3 * pc_y[k] * dsg_85[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pc_y, dsf0_57, dsf0_58, dsf0_59, dsf1_57, \
                         dsf1_58, dsf1_59, dsg_86, dsg_87, dsg_88, \
                         dsg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_11 * dsf0_57[k]
                   - f_12 * dsf1_57[k]
                   + f_3 * pc_y[k] * dsg_86[k];

        t_122[k] = f_6 * dsf0_58[k]
                   - f_7 * dsf1_58[k]
                   + f_3 * pc_y[k] * dsg_87[k];

        t_123[k] = f_4 * dsf0_59[k]
                   - f_5 * dsf1_59[k]
                   + f_3 * pc_y[k] * dsg_88[k];

        t_124[k] = f_3 * pc_y[k] * dsg_89[k];
    }

#pragma omp simd aligned(t_125, pc_z, psg_44, dsf0_59, dsf1_59, \
                         dsg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_0 * psg_44[k]
                   + f_1 * dsf0_59[k]
                   - f_2 * dsf1_59[k]
                   + f_3 * pc_z[k] * dsg_89[k];
    }
}

}  // namespace simdt3ceri
