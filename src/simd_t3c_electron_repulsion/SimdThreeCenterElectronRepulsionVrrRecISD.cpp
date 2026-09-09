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


#include "SimdThreeCenterElectronRepulsionVrrRecISD.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_isd_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsd0,
                                                          const size_t hsp, const size_t hsd1,
                                                          const size_t iss0, const size_t iss1,
                                                          const size_t isp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 2.5 / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 2.0 / q;
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
    auto *t_125 = buffer.data(target + 125);
    auto *t_126 = buffer.data(target + 126);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *hsd0_0 = buffer.data(hsd0 + 0);
    const auto *hsd0_3 = buffer.data(hsd0 + 3);
    const auto *hsd0_5 = buffer.data(hsd0 + 5);
    const auto *hsd0_9 = buffer.data(hsd0 + 9);
    const auto *hsd0_12 = buffer.data(hsd0 + 12);
    const auto *hsd0_17 = buffer.data(hsd0 + 17);
    const auto *hsd0_18 = buffer.data(hsd0 + 18);
    const auto *hsd0_21 = buffer.data(hsd0 + 21);
    const auto *hsd0_30 = buffer.data(hsd0 + 30);
    const auto *hsd0_35 = buffer.data(hsd0 + 35);
    const auto *hsd0_36 = buffer.data(hsd0 + 36);
    const auto *hsd0_39 = buffer.data(hsd0 + 39);
    const auto *hsd0_54 = buffer.data(hsd0 + 54);
    const auto *hsd0_59 = buffer.data(hsd0 + 59);
    const auto *hsd0_60 = buffer.data(hsd0 + 60);
    const auto *hsd0_84 = buffer.data(hsd0 + 84);
    const auto *hsd0_90 = buffer.data(hsd0 + 90);
    const auto *hsd0_93 = buffer.data(hsd0 + 93);
    const auto *hsd0_95 = buffer.data(hsd0 + 95);
    const auto *hsd0_99 = buffer.data(hsd0 + 99);
    const auto *hsd0_101 = buffer.data(hsd0 + 101);
    const auto *hsd0_102 = buffer.data(hsd0 + 102);
    const auto *hsd0_105 = buffer.data(hsd0 + 105);
    const auto *hsd0_107 = buffer.data(hsd0 + 107);
    const auto *hsd0_108 = buffer.data(hsd0 + 108);
    const auto *hsd0_111 = buffer.data(hsd0 + 111);
    const auto *hsd0_113 = buffer.data(hsd0 + 113);
    const auto *hsd0_117 = buffer.data(hsd0 + 117);
    const auto *hsd0_119 = buffer.data(hsd0 + 119);
    const auto *hsd0_120 = buffer.data(hsd0 + 120);
    const auto *hsd0_123 = buffer.data(hsd0 + 123);
    const auto *hsd0_125 = buffer.data(hsd0 + 125);

    const auto *hsp_0 = buffer.data(hsp + 0);
    const auto *hsp_1 = buffer.data(hsp + 1);
    const auto *hsp_2 = buffer.data(hsp + 2);
    const auto *hsp_4 = buffer.data(hsp + 4);
    const auto *hsp_8 = buffer.data(hsp + 8);
    const auto *hsp_9 = buffer.data(hsp + 9);
    const auto *hsp_10 = buffer.data(hsp + 10);
    const auto *hsp_11 = buffer.data(hsp + 11);
    const auto *hsp_13 = buffer.data(hsp + 13);
    const auto *hsp_14 = buffer.data(hsp + 14);
    const auto *hsp_15 = buffer.data(hsp + 15);
    const auto *hsp_16 = buffer.data(hsp + 16);
    const auto *hsp_17 = buffer.data(hsp + 17);
    const auto *hsp_18 = buffer.data(hsp + 18);
    const auto *hsp_19 = buffer.data(hsp + 19);
    const auto *hsp_20 = buffer.data(hsp + 20);
    const auto *hsp_22 = buffer.data(hsp + 22);
    const auto *hsp_23 = buffer.data(hsp + 23);
    const auto *hsp_25 = buffer.data(hsp + 25);
    const auto *hsp_26 = buffer.data(hsp + 26);
    const auto *hsp_27 = buffer.data(hsp + 27);
    const auto *hsp_28 = buffer.data(hsp + 28);
    const auto *hsp_29 = buffer.data(hsp + 29);
    const auto *hsp_30 = buffer.data(hsp + 30);
    const auto *hsp_31 = buffer.data(hsp + 31);
    const auto *hsp_34 = buffer.data(hsp + 34);
    const auto *hsp_35 = buffer.data(hsp + 35);
    const auto *hsp_36 = buffer.data(hsp + 36);
    const auto *hsp_37 = buffer.data(hsp + 37);
    const auto *hsp_38 = buffer.data(hsp + 38);
    const auto *hsp_40 = buffer.data(hsp + 40);
    const auto *hsp_41 = buffer.data(hsp + 41);
    const auto *hsp_42 = buffer.data(hsp + 42);
    const auto *hsp_44 = buffer.data(hsp + 44);
    const auto *hsp_45 = buffer.data(hsp + 45);
    const auto *hsp_46 = buffer.data(hsp + 46);
    const auto *hsp_49 = buffer.data(hsp + 49);
    const auto *hsp_50 = buffer.data(hsp + 50);
    const auto *hsp_51 = buffer.data(hsp + 51);
    const auto *hsp_52 = buffer.data(hsp + 52);
    const auto *hsp_53 = buffer.data(hsp + 53);
    const auto *hsp_54 = buffer.data(hsp + 54);
    const auto *hsp_55 = buffer.data(hsp + 55);
    const auto *hsp_56 = buffer.data(hsp + 56);
    const auto *hsp_58 = buffer.data(hsp + 58);
    const auto *hsp_59 = buffer.data(hsp + 59);
    const auto *hsp_60 = buffer.data(hsp + 60);
    const auto *hsp_62 = buffer.data(hsp + 62);

    const auto *hsd1_0 = buffer.data(hsd1 + 0);
    const auto *hsd1_3 = buffer.data(hsd1 + 3);
    const auto *hsd1_5 = buffer.data(hsd1 + 5);
    const auto *hsd1_9 = buffer.data(hsd1 + 9);
    const auto *hsd1_12 = buffer.data(hsd1 + 12);
    const auto *hsd1_17 = buffer.data(hsd1 + 17);
    const auto *hsd1_18 = buffer.data(hsd1 + 18);
    const auto *hsd1_21 = buffer.data(hsd1 + 21);
    const auto *hsd1_30 = buffer.data(hsd1 + 30);
    const auto *hsd1_35 = buffer.data(hsd1 + 35);
    const auto *hsd1_36 = buffer.data(hsd1 + 36);
    const auto *hsd1_39 = buffer.data(hsd1 + 39);
    const auto *hsd1_54 = buffer.data(hsd1 + 54);
    const auto *hsd1_59 = buffer.data(hsd1 + 59);
    const auto *hsd1_60 = buffer.data(hsd1 + 60);
    const auto *hsd1_84 = buffer.data(hsd1 + 84);
    const auto *hsd1_90 = buffer.data(hsd1 + 90);
    const auto *hsd1_93 = buffer.data(hsd1 + 93);
    const auto *hsd1_95 = buffer.data(hsd1 + 95);
    const auto *hsd1_99 = buffer.data(hsd1 + 99);
    const auto *hsd1_101 = buffer.data(hsd1 + 101);
    const auto *hsd1_102 = buffer.data(hsd1 + 102);
    const auto *hsd1_105 = buffer.data(hsd1 + 105);
    const auto *hsd1_107 = buffer.data(hsd1 + 107);
    const auto *hsd1_108 = buffer.data(hsd1 + 108);
    const auto *hsd1_111 = buffer.data(hsd1 + 111);
    const auto *hsd1_113 = buffer.data(hsd1 + 113);
    const auto *hsd1_117 = buffer.data(hsd1 + 117);
    const auto *hsd1_119 = buffer.data(hsd1 + 119);
    const auto *hsd1_120 = buffer.data(hsd1 + 120);
    const auto *hsd1_123 = buffer.data(hsd1 + 123);
    const auto *hsd1_125 = buffer.data(hsd1 + 125);

    const auto *iss0_0 = buffer.data(iss0 + 0);
    const auto *iss0_1 = buffer.data(iss0 + 1);
    const auto *iss0_2 = buffer.data(iss0 + 2);
    const auto *iss0_3 = buffer.data(iss0 + 3);
    const auto *iss0_5 = buffer.data(iss0 + 5);
    const auto *iss0_6 = buffer.data(iss0 + 6);
    const auto *iss0_7 = buffer.data(iss0 + 7);
    const auto *iss0_8 = buffer.data(iss0 + 8);
    const auto *iss0_9 = buffer.data(iss0 + 9);
    const auto *iss0_10 = buffer.data(iss0 + 10);
    const auto *iss0_11 = buffer.data(iss0 + 11);
    const auto *iss0_12 = buffer.data(iss0 + 12);
    const auto *iss0_13 = buffer.data(iss0 + 13);
    const auto *iss0_14 = buffer.data(iss0 + 14);
    const auto *iss0_21 = buffer.data(iss0 + 21);

    const auto *iss1_0 = buffer.data(iss1 + 0);
    const auto *iss1_1 = buffer.data(iss1 + 1);
    const auto *iss1_2 = buffer.data(iss1 + 2);
    const auto *iss1_3 = buffer.data(iss1 + 3);
    const auto *iss1_5 = buffer.data(iss1 + 5);
    const auto *iss1_6 = buffer.data(iss1 + 6);
    const auto *iss1_7 = buffer.data(iss1 + 7);
    const auto *iss1_8 = buffer.data(iss1 + 8);
    const auto *iss1_9 = buffer.data(iss1 + 9);
    const auto *iss1_10 = buffer.data(iss1 + 10);
    const auto *iss1_11 = buffer.data(iss1 + 11);
    const auto *iss1_12 = buffer.data(iss1 + 12);
    const auto *iss1_13 = buffer.data(iss1 + 13);
    const auto *iss1_14 = buffer.data(iss1 + 14);
    const auto *iss1_21 = buffer.data(iss1 + 21);

    const auto *isp_0 = buffer.data(isp + 0);
    const auto *isp_1 = buffer.data(isp + 1);
    const auto *isp_2 = buffer.data(isp + 2);
    const auto *isp_3 = buffer.data(isp + 3);
    const auto *isp_4 = buffer.data(isp + 4);
    const auto *isp_6 = buffer.data(isp + 6);
    const auto *isp_8 = buffer.data(isp + 8);
    const auto *isp_9 = buffer.data(isp + 9);
    const auto *isp_10 = buffer.data(isp + 10);
    const auto *isp_11 = buffer.data(isp + 11);
    const auto *isp_13 = buffer.data(isp + 13);
    const auto *isp_14 = buffer.data(isp + 14);
    const auto *isp_15 = buffer.data(isp + 15);
    const auto *isp_16 = buffer.data(isp + 16);
    const auto *isp_17 = buffer.data(isp + 17);
    const auto *isp_18 = buffer.data(isp + 18);
    const auto *isp_19 = buffer.data(isp + 19);
    const auto *isp_20 = buffer.data(isp + 20);
    const auto *isp_22 = buffer.data(isp + 22);
    const auto *isp_23 = buffer.data(isp + 23);
    const auto *isp_25 = buffer.data(isp + 25);
    const auto *isp_26 = buffer.data(isp + 26);
    const auto *isp_27 = buffer.data(isp + 27);
    const auto *isp_28 = buffer.data(isp + 28);
    const auto *isp_29 = buffer.data(isp + 29);
    const auto *isp_30 = buffer.data(isp + 30);
    const auto *isp_31 = buffer.data(isp + 31);
    const auto *isp_32 = buffer.data(isp + 32);
    const auto *isp_34 = buffer.data(isp + 34);
    const auto *isp_35 = buffer.data(isp + 35);
    const auto *isp_36 = buffer.data(isp + 36);
    const auto *isp_37 = buffer.data(isp + 37);
    const auto *isp_38 = buffer.data(isp + 38);
    const auto *isp_40 = buffer.data(isp + 40);
    const auto *isp_41 = buffer.data(isp + 41);
    const auto *isp_42 = buffer.data(isp + 42);
    const auto *isp_43 = buffer.data(isp + 43);
    const auto *isp_44 = buffer.data(isp + 44);
    const auto *isp_45 = buffer.data(isp + 45);
    const auto *isp_46 = buffer.data(isp + 46);
    const auto *isp_49 = buffer.data(isp + 49);
    const auto *isp_50 = buffer.data(isp + 50);
    const auto *isp_52 = buffer.data(isp + 52);
    const auto *isp_53 = buffer.data(isp + 53);
    const auto *isp_55 = buffer.data(isp + 55);
    const auto *isp_56 = buffer.data(isp + 56);
    const auto *isp_58 = buffer.data(isp + 58);
    const auto *isp_59 = buffer.data(isp + 59);
    const auto *isp_60 = buffer.data(isp + 60);
    const auto *isp_62 = buffer.data(isp + 62);
    const auto *isp_63 = buffer.data(isp + 63);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, hsp_0, iss0_0, \
                         iss1_0, isp_0, isp_1, isp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hsp_0[k]
                 + f_1 * iss0_0[k]
                 - f_2 * iss1_0[k]
                 + f_3 * pc_x[k] * isp_0[k];

        t_1[k] = f_3 * pc_y[k] * isp_0[k];

        t_2[k] = f_3 * pc_z[k] * isp_0[k];

        t_3[k] = f_1 * iss0_0[k]
                 - f_2 * iss1_0[k]
                 + f_3 * pc_y[k] * isp_1[k];

        t_4[k] = f_3 * pc_y[k] * isp_2[k];

        t_5[k] = f_1 * iss0_0[k]
                 - f_2 * iss1_0[k]
                 + f_3 * pc_z[k] * isp_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pa_y, pc_x, pc_y, pc_z, hsd0_0, hsp_1, hsp_4, \
                         hsd1_0, iss0_1, iss1_1, isp_3, isp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * hsd0_0[k]
                 - f_4 * pc_y[k] * hsd1_0[k];

        t_7[k] = f_5 * hsp_4[k]
                 + f_3 * pc_x[k] * isp_4[k];

        t_8[k] = f_3 * pc_z[k] * isp_3[k];

        t_9[k] = f_6 * hsp_1[k]
                 + f_1 * iss0_1[k]
                 - f_2 * iss1_1[k]
                 + f_3 * pc_y[k] * isp_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pa_y, pa_z, pc_y, pc_z, hsd0_0, hsd0_5, \
                         hsd1_0, hsd1_5, isp_4, isp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_3 * pc_z[k] * isp_4[k];

        t_11[k] = pa_y[k] * hsd0_5[k]
                  - f_4 * pc_y[k] * hsd1_5[k];

        t_12[k] = pa_z[k] * hsd0_0[k]
                  - f_4 * pc_z[k] * hsd1_0[k];

        t_13[k] = f_3 * pc_y[k] * isp_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pa_z, pc_x, pc_y, pc_z, hsd0_3, hsp_2, hsp_8, \
                         hsd1_3, iss0_2, iss1_2, isp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_5 * hsp_8[k]
                  + f_3 * pc_x[k] * isp_8[k];

        t_15[k] = pa_z[k] * hsd0_3[k]
                  - f_4 * pc_z[k] * hsd1_3[k];

        t_16[k] = f_3 * pc_y[k] * isp_8[k];

        t_17[k] = f_6 * hsp_2[k]
                  + f_1 * iss0_2[k]
                  - f_2 * iss1_2[k]
                  + f_3 * pc_z[k] * isp_8[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, pc_x, pc_y, pc_z, hsp_4, hsp_9, hsp_10, \
                         iss0_3, iss1_3, isp_9, isp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_7 * hsp_9[k]
                  + f_1 * iss0_3[k]
                  - f_2 * iss1_3[k]
                  + f_3 * pc_x[k] * isp_9[k];

        t_19[k] = f_7 * hsp_10[k]
                  + f_3 * pc_x[k] * isp_10[k];

        t_20[k] = f_3 * pc_z[k] * isp_9[k];

        t_21[k] = f_8 * hsp_4[k]
                  + f_1 * iss0_3[k]
                  - f_2 * iss1_3[k]
                  + f_3 * pc_y[k] * isp_10[k];

        t_22[k] = f_3 * pc_z[k] * isp_10[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_y, pc_x, pc_y, pc_z, hsd0_12, hsp_13, hsd1_12, \
                         iss0_3, iss1_3, isp_11, isp_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_1 * iss0_3[k]
                  - f_2 * iss1_3[k]
                  + f_3 * pc_z[k] * isp_11[k];

        t_24[k] = pa_y[k] * hsd0_12[k]
                  - f_4 * pc_y[k] * hsd1_12[k];

        t_25[k] = f_7 * hsp_13[k]
                  + f_3 * pc_x[k] * isp_13[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pa_y, pa_z, pc_x, pc_y, pc_z, hsd0_9, \
                         hsd0_17, hsp_8, hsp_14, hsd1_9, hsd1_17, \
                         isp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_7 * hsp_14[k]
                  + f_3 * pc_x[k] * isp_14[k];

        t_27[k] = pa_z[k] * hsd0_9[k]
                  - f_4 * pc_z[k] * hsd1_9[k];

        t_28[k] = f_6 * hsp_8[k]
                  + f_3 * pc_y[k] * isp_14[k];

        t_29[k] = pa_y[k] * hsd0_17[k]
                  - f_4 * pc_y[k] * hsd1_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pc_x, pc_y, hsp_15, hsp_17, iss0_5, \
                         iss1_5, isp_15, isp_16, isp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_7 * hsp_15[k]
                  + f_1 * iss0_5[k]
                  - f_2 * iss1_5[k]
                  + f_3 * pc_x[k] * isp_15[k];

        t_31[k] = f_3 * pc_y[k] * isp_15[k];

        t_32[k] = f_7 * hsp_17[k]
                  + f_3 * pc_x[k] * isp_17[k];

        t_33[k] = f_1 * iss0_5[k]
                  - f_2 * iss1_5[k]
                  + f_3 * pc_y[k] * isp_16[k];

        t_34[k] = f_3 * pc_y[k] * isp_17[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pc_x, pc_z, hsp_8, hsp_18, hsp_19, iss0_5, \
                         iss0_6, iss1_5, iss1_6, isp_17, isp_18, \
                         isp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_8 * hsp_8[k]
                  + f_1 * iss0_5[k]
                  - f_2 * iss1_5[k]
                  + f_3 * pc_z[k] * isp_17[k];

        t_36[k] = f_9 * hsp_18[k]
                  + f_1 * iss0_6[k]
                  - f_2 * iss1_6[k]
                  + f_3 * pc_x[k] * isp_18[k];

        t_37[k] = f_9 * hsp_19[k]
                  + f_3 * pc_x[k] * isp_19[k];

        t_38[k] = f_3 * pc_z[k] * isp_18[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_z, pc_y, pc_z, hsd0_18, hsp_10, hsd1_18, \
                         iss0_6, iss1_6, isp_19, isp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_9 * hsp_10[k]
                  + f_1 * iss0_6[k]
                  - f_2 * iss1_6[k]
                  + f_3 * pc_y[k] * isp_19[k];

        t_40[k] = f_3 * pc_z[k] * isp_19[k];

        t_41[k] = f_1 * iss0_6[k]
                  - f_2 * iss1_6[k]
                  + f_3 * pc_z[k] * isp_20[k];

        t_42[k] = pa_z[k] * hsd0_18[k]
                  - f_4 * pc_z[k] * hsd1_18[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_z, pc_x, pc_y, pc_z, hsd0_21, hsp_14, \
                         hsp_22, hsp_23, hsd1_21, isp_22, isp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_9 * hsp_22[k]
                  + f_3 * pc_x[k] * isp_22[k];

        t_44[k] = f_9 * hsp_23[k]
                  + f_3 * pc_x[k] * isp_23[k];

        t_45[k] = pa_z[k] * hsd0_21[k]
                  - f_4 * pc_z[k] * hsd1_21[k];

        t_46[k] = f_8 * hsp_14[k]
                  + f_3 * pc_y[k] * isp_23[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pa_y, pc_x, pc_y, pc_z, hsd0_30, hsp_11, hsp_25, \
                         hsd1_30, iss0_7, iss1_7, isp_23, isp_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_6 * hsp_11[k]
                  + f_1 * iss0_7[k]
                  - f_2 * iss1_7[k]
                  + f_3 * pc_z[k] * isp_23[k];

        t_48[k] = pa_y[k] * hsd0_30[k]
                  - f_4 * pc_y[k] * hsd1_30[k];

        t_49[k] = f_9 * hsp_25[k]
                  + f_3 * pc_x[k] * isp_25[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_y, pc_x, pc_y, hsd0_35, hsp_16, hsp_17, \
                         hsp_26, hsd1_35, iss0_8, iss1_8, isp_25, \
                         isp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_9 * hsp_26[k]
                  + f_3 * pc_x[k] * isp_26[k];

        t_51[k] = f_6 * hsp_16[k]
                  + f_1 * iss0_8[k]
                  - f_2 * iss1_8[k]
                  + f_3 * pc_y[k] * isp_25[k];

        t_52[k] = f_6 * hsp_17[k]
                  + f_3 * pc_y[k] * isp_26[k];

        t_53[k] = pa_y[k] * hsd0_35[k]
                  - f_4 * pc_y[k] * hsd1_35[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pc_x, pc_y, hsp_27, hsp_29, iss0_9, \
                         iss1_9, isp_27, isp_28, isp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_9 * hsp_27[k]
                  + f_1 * iss0_9[k]
                  - f_2 * iss1_9[k]
                  + f_3 * pc_x[k] * isp_27[k];

        t_55[k] = f_3 * pc_y[k] * isp_27[k];

        t_56[k] = f_9 * hsp_29[k]
                  + f_3 * pc_x[k] * isp_29[k];

        t_57[k] = f_1 * iss0_9[k]
                  - f_2 * iss1_9[k]
                  + f_3 * pc_y[k] * isp_28[k];

        t_58[k] = f_3 * pc_y[k] * isp_29[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pc_x, pc_z, hsp_17, hsp_30, hsp_31, iss0_9, \
                         iss0_10, iss1_9, iss1_10, isp_29, isp_30, \
                         isp_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_9 * hsp_17[k]
                  + f_1 * iss0_9[k]
                  - f_2 * iss1_9[k]
                  + f_3 * pc_z[k] * isp_29[k];

        t_60[k] = f_8 * hsp_30[k]
                  + f_1 * iss0_10[k]
                  - f_2 * iss1_10[k]
                  + f_3 * pc_x[k] * isp_30[k];

        t_61[k] = f_8 * hsp_31[k]
                  + f_3 * pc_x[k] * isp_31[k];

        t_62[k] = f_3 * pc_z[k] * isp_30[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, pa_z, pc_y, pc_z, hsd0_36, hsp_19, hsd1_36, \
                         iss0_10, iss1_10, isp_31, isp_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_7 * hsp_19[k]
                  + f_1 * iss0_10[k]
                  - f_2 * iss1_10[k]
                  + f_3 * pc_y[k] * isp_31[k];

        t_64[k] = f_3 * pc_z[k] * isp_31[k];

        t_65[k] = f_1 * iss0_10[k]
                  - f_2 * iss1_10[k]
                  + f_3 * pc_z[k] * isp_32[k];

        t_66[k] = pa_z[k] * hsd0_36[k]
                  - f_4 * pc_z[k] * hsd1_36[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, t_70, pa_z, pc_x, pc_y, pc_z, hsd0_39, hsp_23, \
                         hsp_34, hsp_35, hsd1_39, isp_34, isp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_8 * hsp_34[k]
                  + f_3 * pc_x[k] * isp_34[k];

        t_68[k] = f_8 * hsp_35[k]
                  + f_3 * pc_x[k] * isp_35[k];

        t_69[k] = pa_z[k] * hsd0_39[k]
                  - f_4 * pc_z[k] * hsd1_39[k];

        t_70[k] = f_9 * hsp_23[k]
                  + f_3 * pc_y[k] * isp_35[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, pc_x, pc_z, hsp_20, hsp_36, hsp_37, iss0_11, \
                         iss0_12, iss1_11, iss1_12, isp_35, isp_36, \
                         isp_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_6 * hsp_20[k]
                  + f_1 * iss0_11[k]
                  - f_2 * iss1_11[k]
                  + f_3 * pc_z[k] * isp_35[k];

        t_72[k] = f_8 * hsp_36[k]
                  + f_1 * iss0_12[k]
                  - f_2 * iss1_12[k]
                  + f_3 * pc_x[k] * isp_36[k];

        t_73[k] = f_8 * hsp_37[k]
                  + f_3 * pc_x[k] * isp_37[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pc_x, pc_y, pc_z, hsp_23, hsp_25, hsp_26, \
                         hsp_38, iss0_12, iss1_12, isp_37, isp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_8 * hsp_38[k]
                  + f_3 * pc_x[k] * isp_38[k];

        t_75[k] = f_8 * hsp_25[k]
                  + f_1 * iss0_12[k]
                  - f_2 * iss1_12[k]
                  + f_3 * pc_y[k] * isp_37[k];

        t_76[k] = f_8 * hsp_26[k]
                  + f_3 * pc_y[k] * isp_38[k];

        t_77[k] = f_8 * hsp_23[k]
                  + f_1 * iss0_12[k]
                  - f_2 * iss1_12[k]
                  + f_3 * pc_z[k] * isp_38[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_y, pc_x, pc_y, hsd0_54, hsp_28, hsp_40, \
                         hsp_41, hsd1_54, iss0_13, iss1_13, isp_40, \
                         isp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = pa_y[k] * hsd0_54[k]
                  - f_4 * pc_y[k] * hsd1_54[k];

        t_79[k] = f_8 * hsp_40[k]
                  + f_3 * pc_x[k] * isp_40[k];

        t_80[k] = f_8 * hsp_41[k]
                  + f_3 * pc_x[k] * isp_41[k];

        t_81[k] = f_6 * hsp_28[k]
                  + f_1 * iss0_13[k]
                  - f_2 * iss1_13[k]
                  + f_3 * pc_y[k] * isp_40[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pa_y, pc_x, pc_y, hsd0_59, hsp_29, hsp_42, \
                         hsd1_59, iss0_14, iss1_14, isp_41, isp_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_6 * hsp_29[k]
                  + f_3 * pc_y[k] * isp_41[k];

        t_83[k] = pa_y[k] * hsd0_59[k]
                  - f_4 * pc_y[k] * hsd1_59[k];

        t_84[k] = f_8 * hsp_42[k]
                  + f_1 * iss0_14[k]
                  - f_2 * iss1_14[k]
                  + f_3 * pc_x[k] * isp_42[k];

        t_85[k] = f_3 * pc_y[k] * isp_42[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pc_x, pc_y, pc_z, hsp_29, hsp_44, iss0_14, \
                         iss1_14, isp_43, isp_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_8 * hsp_44[k]
                  + f_3 * pc_x[k] * isp_44[k];

        t_87[k] = f_1 * iss0_14[k]
                  - f_2 * iss1_14[k]
                  + f_3 * pc_y[k] * isp_43[k];

        t_88[k] = f_3 * pc_y[k] * isp_44[k];

        t_89[k] = f_7 * hsp_29[k]
                  + f_1 * iss0_14[k]
                  - f_2 * iss1_14[k]
                  + f_3 * pc_z[k] * isp_44[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, pa_x, pc_x, pc_z, hsd0_90, hsd0_93, \
                         hsp_45, hsp_46, hsd1_90, hsd1_93, isp_45, \
                         isp_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = pa_x[k] * hsd0_90[k]
                  + f_8 * hsp_45[k]
                  - f_4 * pc_x[k] * hsd1_90[k];

        t_91[k] = f_6 * hsp_46[k]
                  + f_3 * pc_x[k] * isp_46[k];

        t_92[k] = f_3 * pc_z[k] * isp_45[k];

        t_93[k] = pa_x[k] * hsd0_93[k]
                  - f_4 * pc_x[k] * hsd1_93[k];

        t_94[k] = f_3 * pc_z[k] * isp_46[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pa_x, pa_z, pc_x, pc_z, hsd0_60, hsd0_95, \
                         hsp_49, hsp_50, hsd1_60, hsd1_95, isp_49, \
                         isp_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = pa_x[k] * hsd0_95[k]
                  - f_4 * pc_x[k] * hsd1_95[k];

        t_96[k] = pa_z[k] * hsd0_60[k]
                  - f_4 * pc_z[k] * hsd1_60[k];

        t_97[k] = f_6 * hsp_49[k]
                  + f_3 * pc_x[k] * isp_49[k];

        t_98[k] = f_6 * hsp_50[k]
                  + f_3 * pc_x[k] * isp_50[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pa_x, pc_x, pc_y, hsd0_99, hsd0_101, \
                         hsd0_102, hsp_35, hsp_51, hsd1_99, hsd1_101, hsd1_102, \
                         isp_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = pa_x[k] * hsd0_99[k]
                  - f_4 * pc_x[k] * hsd1_99[k];

        t_100[k] = f_7 * hsp_35[k]
                   + f_3 * pc_y[k] * isp_50[k];

        t_101[k] = pa_x[k] * hsd0_101[k]
                   - f_4 * pc_x[k] * hsd1_101[k];

        t_102[k] = pa_x[k] * hsd0_102[k]
                   + f_8 * hsp_51[k]
                   - f_4 * pc_x[k] * hsd1_102[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, pa_x, pc_x, pc_y, hsd0_105, hsp_38, \
                         hsp_52, hsp_53, hsd1_105, isp_52, isp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_6 * hsp_52[k]
                   + f_3 * pc_x[k] * isp_52[k];

        t_104[k] = f_6 * hsp_53[k]
                   + f_3 * pc_x[k] * isp_53[k];

        t_105[k] = pa_x[k] * hsd0_105[k]
                   - f_4 * pc_x[k] * hsd1_105[k];

        t_106[k] = f_9 * hsp_38[k]
                   + f_3 * pc_y[k] * isp_53[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pa_x, pc_x, hsd0_107, hsd0_108, hsp_54, \
                         hsp_55, hsp_56, hsd1_107, hsd1_108, isp_55, \
                         isp_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = pa_x[k] * hsd0_107[k]
                   - f_4 * pc_x[k] * hsd1_107[k];

        t_108[k] = pa_x[k] * hsd0_108[k]
                   + f_8 * hsp_54[k]
                   - f_4 * pc_x[k] * hsd1_108[k];

        t_109[k] = f_6 * hsp_55[k]
                   + f_3 * pc_x[k] * isp_55[k];

        t_110[k] = f_6 * hsp_56[k]
                   + f_3 * pc_x[k] * isp_56[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pa_x, pa_y, pc_x, pc_y, hsd0_84, \
                         hsd0_111, hsd0_113, hsp_41, hsd1_84, hsd1_111, hsd1_113, \
                         isp_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = pa_x[k] * hsd0_111[k]
                   - f_4 * pc_x[k] * hsd1_111[k];

        t_112[k] = f_8 * hsp_41[k]
                   + f_3 * pc_y[k] * isp_56[k];

        t_113[k] = pa_x[k] * hsd0_113[k]
                   - f_4 * pc_x[k] * hsd1_113[k];

        t_114[k] = pa_y[k] * hsd0_84[k]
                   - f_4 * pc_y[k] * hsd1_84[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, pa_x, pc_x, pc_y, hsd0_117, hsp_44, \
                         hsp_58, hsp_59, hsd1_117, isp_58, isp_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_6 * hsp_58[k]
                   + f_3 * pc_x[k] * isp_58[k];

        t_116[k] = f_6 * hsp_59[k]
                   + f_3 * pc_x[k] * isp_59[k];

        t_117[k] = pa_x[k] * hsd0_117[k]
                   - f_4 * pc_x[k] * hsd1_117[k];

        t_118[k] = f_6 * hsp_44[k]
                   + f_3 * pc_y[k] * isp_59[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, pa_x, pc_x, pc_y, hsd0_119, hsd0_120, \
                         hsp_60, hsp_62, hsd1_119, hsd1_120, isp_60, \
                         isp_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = pa_x[k] * hsd0_119[k]
                   - f_4 * pc_x[k] * hsd1_119[k];

        t_120[k] = pa_x[k] * hsd0_120[k]
                   + f_8 * hsp_60[k]
                   - f_4 * pc_x[k] * hsd1_120[k];

        t_121[k] = f_3 * pc_y[k] * isp_60[k];

        t_122[k] = f_6 * hsp_62[k]
                   + f_3 * pc_x[k] * isp_62[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, pa_x, pc_x, pc_y, hsd0_123, hsd0_125, \
                         hsd1_123, hsd1_125, iss0_21, iss1_21, isp_62, \
                         isp_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = pa_x[k] * hsd0_123[k]
                   - f_4 * pc_x[k] * hsd1_123[k];

        t_124[k] = f_3 * pc_y[k] * isp_62[k];

        t_125[k] = pa_x[k] * hsd0_125[k]
                   - f_4 * pc_x[k] * hsd1_125[k];

        t_126[k] = f_1 * iss0_21[k]
                   - f_2 * iss1_21[k]
                   + f_3 * pc_x[k] * isp_63[k];
    }
}

static auto
compute_prim_isd_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t hsd0,
                                                          const size_t hsp, const size_t hsd1,
                                                          const size_t iss0, const size_t iss1,
                                                          const size_t isp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 2.5 / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 2.0 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 1.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *hsd0_90 = buffer.data(hsd0 + 90);
    const auto *hsd0_93 = buffer.data(hsd0 + 93);
    const auto *hsd0_120 = buffer.data(hsd0 + 120);
    const auto *hsd0_123 = buffer.data(hsd0 + 123);
    const auto *hsd0_125 = buffer.data(hsd0 + 125);

    const auto *hsp_46 = buffer.data(hsp + 46);
    const auto *hsp_47 = buffer.data(hsp + 47);
    const auto *hsp_50 = buffer.data(hsp + 50);
    const auto *hsp_52 = buffer.data(hsp + 52);
    const auto *hsp_53 = buffer.data(hsp + 53);
    const auto *hsp_55 = buffer.data(hsp + 55);
    const auto *hsp_56 = buffer.data(hsp + 56);
    const auto *hsp_58 = buffer.data(hsp + 58);
    const auto *hsp_59 = buffer.data(hsp + 59);
    const auto *hsp_61 = buffer.data(hsp + 61);
    const auto *hsp_62 = buffer.data(hsp + 62);

    const auto *hsd1_90 = buffer.data(hsd1 + 90);
    const auto *hsd1_93 = buffer.data(hsd1 + 93);
    const auto *hsd1_120 = buffer.data(hsd1 + 120);
    const auto *hsd1_123 = buffer.data(hsd1 + 123);
    const auto *hsd1_125 = buffer.data(hsd1 + 125);

    const auto *iss0_21 = buffer.data(iss0 + 21);
    const auto *iss0_22 = buffer.data(iss0 + 22);
    const auto *iss0_23 = buffer.data(iss0 + 23);
    const auto *iss0_24 = buffer.data(iss0 + 24);
    const auto *iss0_25 = buffer.data(iss0 + 25);
    const auto *iss0_27 = buffer.data(iss0 + 27);

    const auto *iss1_21 = buffer.data(iss1 + 21);
    const auto *iss1_22 = buffer.data(iss1 + 22);
    const auto *iss1_23 = buffer.data(iss1 + 23);
    const auto *iss1_24 = buffer.data(iss1 + 24);
    const auto *iss1_25 = buffer.data(iss1 + 25);
    const auto *iss1_27 = buffer.data(iss1 + 27);

    const auto *isp_64 = buffer.data(isp + 64);
    const auto *isp_65 = buffer.data(isp + 65);
    const auto *isp_67 = buffer.data(isp + 67);
    const auto *isp_68 = buffer.data(isp + 68);
    const auto *isp_69 = buffer.data(isp + 69);
    const auto *isp_70 = buffer.data(isp + 70);
    const auto *isp_71 = buffer.data(isp + 71);
    const auto *isp_72 = buffer.data(isp + 72);
    const auto *isp_73 = buffer.data(isp + 73);
    const auto *isp_74 = buffer.data(isp + 74);
    const auto *isp_75 = buffer.data(isp + 75);
    const auto *isp_76 = buffer.data(isp + 76);
    const auto *isp_77 = buffer.data(isp + 77);
    const auto *isp_79 = buffer.data(isp + 79);
    const auto *isp_80 = buffer.data(isp + 80);
    const auto *isp_81 = buffer.data(isp + 81);
    const auto *isp_82 = buffer.data(isp + 82);
    const auto *isp_83 = buffer.data(isp + 83);

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, pc_x, pc_y, pc_z, hsp_46, iss0_21, \
                         iss1_21, isp_64, isp_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_3 * pc_x[k] * isp_64[k];

        t_128[k] = f_3 * pc_x[k] * isp_65[k];

        t_129[k] = f_0 * hsp_46[k]
                   + f_1 * iss0_21[k]
                   - f_2 * iss1_21[k]
                   + f_3 * pc_y[k] * isp_64[k];

        t_130[k] = f_3 * pc_z[k] * isp_64[k];

        t_131[k] = f_1 * iss0_21[k]
                   - f_2 * iss1_21[k]
                   + f_3 * pc_z[k] * isp_65[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, pa_z, pc_x, pc_y, pc_z, hsd0_90, \
                         hsd0_93, hsp_50, hsd1_90, hsd1_93, isp_67, \
                         isp_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = pa_z[k] * hsd0_90[k]
                   - f_4 * pc_z[k] * hsd1_90[k];

        t_133[k] = f_3 * pc_x[k] * isp_67[k];

        t_134[k] = f_3 * pc_x[k] * isp_68[k];

        t_135[k] = pa_z[k] * hsd0_93[k]
                   - f_4 * pc_z[k] * hsd1_93[k];

        t_136[k] = f_5 * hsp_50[k]
                   + f_3 * pc_y[k] * isp_68[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, pc_x, pc_z, hsp_47, iss0_22, iss0_23, \
                         iss1_22, iss1_23, isp_68, isp_69, isp_70, \
                         isp_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_6 * hsp_47[k]
                   + f_1 * iss0_22[k]
                   - f_2 * iss1_22[k]
                   + f_3 * pc_z[k] * isp_68[k];

        t_138[k] = f_1 * iss0_23[k]
                   - f_2 * iss1_23[k]
                   + f_3 * pc_x[k] * isp_69[k];

        t_139[k] = f_3 * pc_x[k] * isp_70[k];

        t_140[k] = f_3 * pc_x[k] * isp_71[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, pc_y, pc_z, hsp_50, hsp_52, hsp_53, iss0_23, \
                         iss1_23, isp_70, isp_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_7 * hsp_52[k]
                   + f_1 * iss0_23[k]
                   - f_2 * iss1_23[k]
                   + f_3 * pc_y[k] * isp_70[k];

        t_142[k] = f_7 * hsp_53[k]
                   + f_3 * pc_y[k] * isp_71[k];

        t_143[k] = f_8 * hsp_50[k]
                   + f_1 * iss0_23[k]
                   - f_2 * iss1_23[k]
                   + f_3 * pc_z[k] * isp_71[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, pc_x, pc_y, hsp_55, hsp_56, \
                         iss0_24, iss1_24, isp_72, isp_73, isp_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_1 * iss0_24[k]
                   - f_2 * iss1_24[k]
                   + f_3 * pc_x[k] * isp_72[k];

        t_145[k] = f_3 * pc_x[k] * isp_73[k];

        t_146[k] = f_3 * pc_x[k] * isp_74[k];

        t_147[k] = f_9 * hsp_55[k]
                   + f_1 * iss0_24[k]
                   - f_2 * iss1_24[k]
                   + f_3 * pc_y[k] * isp_73[k];

        t_148[k] = f_9 * hsp_56[k]
                   + f_3 * pc_y[k] * isp_74[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, pc_x, pc_z, hsp_53, iss0_24, iss0_25, \
                         iss1_24, iss1_25, isp_74, isp_75, isp_76, \
                         isp_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_9 * hsp_53[k]
                   + f_1 * iss0_24[k]
                   - f_2 * iss1_24[k]
                   + f_3 * pc_z[k] * isp_74[k];

        t_150[k] = f_1 * iss0_25[k]
                   - f_2 * iss1_25[k]
                   + f_3 * pc_x[k] * isp_75[k];

        t_151[k] = f_3 * pc_x[k] * isp_76[k];

        t_152[k] = f_3 * pc_x[k] * isp_77[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, pa_y, pc_y, pc_z, hsd0_120, hsp_56, \
                         hsp_58, hsp_59, hsd1_120, iss0_25, iss1_25, isp_76, \
                         isp_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_8 * hsp_58[k]
                   + f_1 * iss0_25[k]
                   - f_2 * iss1_25[k]
                   + f_3 * pc_y[k] * isp_76[k];

        t_154[k] = f_8 * hsp_59[k]
                   + f_3 * pc_y[k] * isp_77[k];

        t_155[k] = f_7 * hsp_56[k]
                   + f_1 * iss0_25[k]
                   - f_2 * iss1_25[k]
                   + f_3 * pc_z[k] * isp_77[k];

        t_156[k] = pa_y[k] * hsd0_120[k]
                   - f_4 * pc_y[k] * hsd1_120[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, pa_y, pc_x, pc_y, hsd0_123, \
                         hsd0_125, hsp_61, hsp_62, hsd1_123, hsd1_125, isp_79, \
                         isp_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_3 * pc_x[k] * isp_79[k];

        t_158[k] = f_3 * pc_x[k] * isp_80[k];

        t_159[k] = pa_y[k] * hsd0_123[k]
                   + f_8 * hsp_61[k]
                   - f_4 * pc_y[k] * hsd1_123[k];

        t_160[k] = f_6 * hsp_62[k]
                   + f_3 * pc_y[k] * isp_80[k];

        t_161[k] = pa_y[k] * hsd0_125[k]
                   - f_4 * pc_y[k] * hsd1_125[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, t_167, pc_x, pc_y, pc_z, hsp_62, \
                         iss0_27, iss1_27, isp_81, isp_82, isp_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_1 * iss0_27[k]
                   - f_2 * iss1_27[k]
                   + f_3 * pc_x[k] * isp_81[k];

        t_163[k] = f_3 * pc_x[k] * isp_82[k];

        t_164[k] = f_3 * pc_x[k] * isp_83[k];

        t_165[k] = f_1 * iss0_27[k]
                   - f_2 * iss1_27[k]
                   + f_3 * pc_y[k] * isp_82[k];

        t_166[k] = f_3 * pc_y[k] * isp_83[k];

        t_167[k] = f_0 * hsp_62[k]
                   + f_1 * iss0_27[k]
                   - f_2 * iss1_27[k]
                   + f_3 * pc_z[k] * isp_83[k];
    }
}

auto
compute_prim_isd_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t hsd0, const size_t hsp,
                                                   const size_t hsd1, const size_t iss0,
                                                   const size_t iss1, const size_t isp,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_isd_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, hsd0, hsp,
                                                              hsd1, iss0, iss1, isp, ncols,
                                                              gamma, p, q);

    compute_prim_isd_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, hsd0, hsp,
                                                              hsd1, iss0, iss1, isp, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
