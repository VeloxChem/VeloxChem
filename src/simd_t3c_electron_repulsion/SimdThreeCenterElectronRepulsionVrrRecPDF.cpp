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


#include "SimdThreeCenterElectronRepulsionVrrRecPDF.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_pdf_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sdf0, const size_t sdd,
                                                          const size_t sdf1, const size_t ppf0,
                                                          const size_t ppd, const size_t ppf1,
                                                          const size_t pdp0, const size_t pdp1,
                                                          const size_t pdd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.0 / q;
    const auto f_2 = 1.0 / gamma;
    const auto f_3 = p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = gamma / q;
    const auto f_6 = 1.5 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sdf0_0 = buffer.data(sdf0 + 0);
    const auto *sdf0_1 = buffer.data(sdf0 + 1);
    const auto *sdf0_2 = buffer.data(sdf0 + 2);
    const auto *sdf0_6 = buffer.data(sdf0 + 6);
    const auto *sdf0_9 = buffer.data(sdf0 + 9);
    const auto *sdf0_20 = buffer.data(sdf0 + 20);
    const auto *sdf0_29 = buffer.data(sdf0 + 29);
    const auto *sdf0_30 = buffer.data(sdf0 + 30);
    const auto *sdf0_36 = buffer.data(sdf0 + 36);
    const auto *sdf0_39 = buffer.data(sdf0 + 39);
    const auto *sdf0_46 = buffer.data(sdf0 + 46);
    const auto *sdf0_49 = buffer.data(sdf0 + 49);
    const auto *sdf0_50 = buffer.data(sdf0 + 50);
    const auto *sdf0_56 = buffer.data(sdf0 + 56);
    const auto *sdf0_59 = buffer.data(sdf0 + 59);

    const auto *sdd_0 = buffer.data(sdd + 0);
    const auto *sdd_3 = buffer.data(sdd + 3);
    const auto *sdd_5 = buffer.data(sdd + 5);
    const auto *sdd_9 = buffer.data(sdd + 9);
    const auto *sdd_17 = buffer.data(sdd + 17);
    const auto *sdd_18 = buffer.data(sdd + 18);
    const auto *sdd_21 = buffer.data(sdd + 21);
    const auto *sdd_23 = buffer.data(sdd + 23);
    const auto *sdd_27 = buffer.data(sdd + 27);
    const auto *sdd_29 = buffer.data(sdd + 29);
    const auto *sdd_30 = buffer.data(sdd + 30);
    const auto *sdd_33 = buffer.data(sdd + 33);
    const auto *sdd_35 = buffer.data(sdd + 35);

    const auto *sdf1_0 = buffer.data(sdf1 + 0);
    const auto *sdf1_1 = buffer.data(sdf1 + 1);
    const auto *sdf1_2 = buffer.data(sdf1 + 2);
    const auto *sdf1_6 = buffer.data(sdf1 + 6);
    const auto *sdf1_9 = buffer.data(sdf1 + 9);
    const auto *sdf1_20 = buffer.data(sdf1 + 20);
    const auto *sdf1_29 = buffer.data(sdf1 + 29);
    const auto *sdf1_30 = buffer.data(sdf1 + 30);
    const auto *sdf1_36 = buffer.data(sdf1 + 36);
    const auto *sdf1_39 = buffer.data(sdf1 + 39);
    const auto *sdf1_46 = buffer.data(sdf1 + 46);
    const auto *sdf1_49 = buffer.data(sdf1 + 49);
    const auto *sdf1_50 = buffer.data(sdf1 + 50);
    const auto *sdf1_56 = buffer.data(sdf1 + 56);
    const auto *sdf1_59 = buffer.data(sdf1 + 59);

    const auto *ppf0_0 = buffer.data(ppf0 + 0);
    const auto *ppf0_3 = buffer.data(ppf0 + 3);
    const auto *ppf0_5 = buffer.data(ppf0 + 5);
    const auto *ppf0_20 = buffer.data(ppf0 + 20);
    const auto *ppf0_31 = buffer.data(ppf0 + 31);
    const auto *ppf0_41 = buffer.data(ppf0 + 41);
    const auto *ppf0_46 = buffer.data(ppf0 + 46);
    const auto *ppf0_48 = buffer.data(ppf0 + 48);
    const auto *ppf0_49 = buffer.data(ppf0 + 49);
    const auto *ppf0_56 = buffer.data(ppf0 + 56);

    const auto *ppd_0 = buffer.data(ppd + 0);
    const auto *ppd_2 = buffer.data(ppd + 2);
    const auto *ppd_3 = buffer.data(ppd + 3);
    const auto *ppd_5 = buffer.data(ppd + 5);
    const auto *ppd_6 = buffer.data(ppd + 6);
    const auto *ppd_8 = buffer.data(ppd + 8);
    const auto *ppd_9 = buffer.data(ppd + 9);
    const auto *ppd_11 = buffer.data(ppd + 11);
    const auto *ppd_12 = buffer.data(ppd + 12);
    const auto *ppd_14 = buffer.data(ppd + 14);
    const auto *ppd_15 = buffer.data(ppd + 15);
    const auto *ppd_17 = buffer.data(ppd + 17);
    const auto *ppd_18 = buffer.data(ppd + 18);
    const auto *ppd_21 = buffer.data(ppd + 21);
    const auto *ppd_22 = buffer.data(ppd + 22);
    const auto *ppd_23 = buffer.data(ppd + 23);
    const auto *ppd_24 = buffer.data(ppd + 24);
    const auto *ppd_25 = buffer.data(ppd + 25);
    const auto *ppd_27 = buffer.data(ppd + 27);
    const auto *ppd_28 = buffer.data(ppd + 28);
    const auto *ppd_29 = buffer.data(ppd + 29);
    const auto *ppd_30 = buffer.data(ppd + 30);
    const auto *ppd_33 = buffer.data(ppd + 33);
    const auto *ppd_34 = buffer.data(ppd + 34);
    const auto *ppd_35 = buffer.data(ppd + 35);
    const auto *ppd_39 = buffer.data(ppd + 39);
    const auto *ppd_40 = buffer.data(ppd + 40);
    const auto *ppd_41 = buffer.data(ppd + 41);

    const auto *ppf1_0 = buffer.data(ppf1 + 0);
    const auto *ppf1_3 = buffer.data(ppf1 + 3);
    const auto *ppf1_5 = buffer.data(ppf1 + 5);
    const auto *ppf1_20 = buffer.data(ppf1 + 20);
    const auto *ppf1_31 = buffer.data(ppf1 + 31);
    const auto *ppf1_41 = buffer.data(ppf1 + 41);
    const auto *ppf1_46 = buffer.data(ppf1 + 46);
    const auto *ppf1_48 = buffer.data(ppf1 + 48);
    const auto *ppf1_49 = buffer.data(ppf1 + 49);
    const auto *ppf1_56 = buffer.data(ppf1 + 56);

    const auto *pdp0_0 = buffer.data(pdp0 + 0);
    const auto *pdp0_1 = buffer.data(pdp0 + 1);
    const auto *pdp0_2 = buffer.data(pdp0 + 2);
    const auto *pdp0_4 = buffer.data(pdp0 + 4);
    const auto *pdp0_5 = buffer.data(pdp0 + 5);
    const auto *pdp0_7 = buffer.data(pdp0 + 7);
    const auto *pdp0_8 = buffer.data(pdp0 + 8);
    const auto *pdp0_21 = buffer.data(pdp0 + 21);
    const auto *pdp0_27 = buffer.data(pdp0 + 27);
    const auto *pdp0_28 = buffer.data(pdp0 + 28);
    const auto *pdp0_29 = buffer.data(pdp0 + 29);
    const auto *pdp0_30 = buffer.data(pdp0 + 30);
    const auto *pdp0_32 = buffer.data(pdp0 + 32);
    const auto *pdp0_34 = buffer.data(pdp0 + 34);
    const auto *pdp0_38 = buffer.data(pdp0 + 38);

    const auto *pdp1_0 = buffer.data(pdp1 + 0);
    const auto *pdp1_1 = buffer.data(pdp1 + 1);
    const auto *pdp1_2 = buffer.data(pdp1 + 2);
    const auto *pdp1_4 = buffer.data(pdp1 + 4);
    const auto *pdp1_5 = buffer.data(pdp1 + 5);
    const auto *pdp1_7 = buffer.data(pdp1 + 7);
    const auto *pdp1_8 = buffer.data(pdp1 + 8);
    const auto *pdp1_21 = buffer.data(pdp1 + 21);
    const auto *pdp1_27 = buffer.data(pdp1 + 27);
    const auto *pdp1_28 = buffer.data(pdp1 + 28);
    const auto *pdp1_29 = buffer.data(pdp1 + 29);
    const auto *pdp1_30 = buffer.data(pdp1 + 30);
    const auto *pdp1_32 = buffer.data(pdp1 + 32);
    const auto *pdp1_34 = buffer.data(pdp1 + 34);
    const auto *pdp1_38 = buffer.data(pdp1 + 38);

    const auto *pdd_0 = buffer.data(pdd + 0);
    const auto *pdd_2 = buffer.data(pdd + 2);
    const auto *pdd_3 = buffer.data(pdd + 3);
    const auto *pdd_5 = buffer.data(pdd + 5);
    const auto *pdd_6 = buffer.data(pdd + 6);
    const auto *pdd_8 = buffer.data(pdd + 8);
    const auto *pdd_9 = buffer.data(pdd + 9);
    const auto *pdd_11 = buffer.data(pdd + 11);
    const auto *pdd_12 = buffer.data(pdd + 12);
    const auto *pdd_14 = buffer.data(pdd + 14);
    const auto *pdd_15 = buffer.data(pdd + 15);
    const auto *pdd_17 = buffer.data(pdd + 17);
    const auto *pdd_18 = buffer.data(pdd + 18);
    const auto *pdd_20 = buffer.data(pdd + 20);
    const auto *pdd_21 = buffer.data(pdd + 21);
    const auto *pdd_23 = buffer.data(pdd + 23);
    const auto *pdd_24 = buffer.data(pdd + 24);
    const auto *pdd_26 = buffer.data(pdd + 26);
    const auto *pdd_27 = buffer.data(pdd + 27);
    const auto *pdd_29 = buffer.data(pdd + 29);
    const auto *pdd_30 = buffer.data(pdd + 30);
    const auto *pdd_32 = buffer.data(pdd + 32);
    const auto *pdd_33 = buffer.data(pdd + 33);
    const auto *pdd_35 = buffer.data(pdd + 35);
    const auto *pdd_36 = buffer.data(pdd + 36);
    const auto *pdd_39 = buffer.data(pdd + 39);
    const auto *pdd_40 = buffer.data(pdd + 40);
    const auto *pdd_41 = buffer.data(pdd + 41);
    const auto *pdd_42 = buffer.data(pdd + 42);
    const auto *pdd_45 = buffer.data(pdd + 45);
    const auto *pdd_46 = buffer.data(pdd + 46);
    const auto *pdd_47 = buffer.data(pdd + 47);
    const auto *pdd_48 = buffer.data(pdd + 48);
    const auto *pdd_51 = buffer.data(pdd + 51);
    const auto *pdd_52 = buffer.data(pdd + 52);
    const auto *pdd_53 = buffer.data(pdd + 53);
    const auto *pdd_54 = buffer.data(pdd + 54);
    const auto *pdd_55 = buffer.data(pdd + 55);
    const auto *pdd_57 = buffer.data(pdd + 57);
    const auto *pdd_58 = buffer.data(pdd + 58);
    const auto *pdd_59 = buffer.data(pdd + 59);
    const auto *pdd_60 = buffer.data(pdd + 60);
    const auto *pdd_63 = buffer.data(pdd + 63);
    const auto *pdd_64 = buffer.data(pdd + 64);
    const auto *pdd_65 = buffer.data(pdd + 65);
    const auto *pdd_66 = buffer.data(pdd + 66);
    const auto *pdd_67 = buffer.data(pdd + 67);
    const auto *pdd_69 = buffer.data(pdd + 69);
    const auto *pdd_70 = buffer.data(pdd + 70);
    const auto *pdd_71 = buffer.data(pdd + 71);
    const auto *pdd_72 = buffer.data(pdd + 72);
    const auto *pdd_75 = buffer.data(pdd + 75);
    const auto *pdd_76 = buffer.data(pdd + 76);
    const auto *pdd_77 = buffer.data(pdd + 77);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, sdd_0, sdd_3, ppd_0, ppd_3, \
                         pdp0_0, pdp1_0, pdd_0, pdd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sdd_0[k]
                 + f_1 * ppd_0[k]
                 + f_2 * pdp0_0[k]
                 - f_3 * pdp1_0[k]
                 + f_4 * pc_x[k] * pdd_0[k];

        t_1[k] = f_4 * pc_y[k] * pdd_0[k];

        t_2[k] = f_4 * pc_z[k] * pdd_0[k];

        t_3[k] = f_0 * sdd_3[k]
                 + f_1 * ppd_3[k]
                 + f_4 * pc_x[k] * pdd_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pc_x, pc_y, pc_z, sdd_5, ppd_5, pdp0_1, \
                         pdp1_1, pdd_2, pdd_3, pdd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_4 * pc_y[k] * pdd_2[k];

        t_5[k] = f_0 * sdd_5[k]
                 + f_1 * ppd_5[k]
                 + f_4 * pc_x[k] * pdd_5[k];

        t_6[k] = f_2 * pdp0_1[k]
                 - f_3 * pdp1_1[k]
                 + f_4 * pc_y[k] * pdd_3[k];

        t_7[k] = f_4 * pc_z[k] * pdd_3[k];

        t_8[k] = f_4 * pc_y[k] * pdd_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pb_y, pc_y, pc_z, ppf0_0, ppd_0, ppf1_0, \
                         pdp0_2, pdp1_2, pdd_5, pdd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_2 * pdp0_2[k]
                 - f_3 * pdp1_2[k]
                 + f_4 * pc_z[k] * pdd_5[k];

        t_10[k] = pb_y[k] * ppf0_0[k]
                  - f_5 * pc_y[k] * ppf1_0[k];

        t_11[k] = f_0 * ppd_0[k]
                  + f_4 * pc_y[k] * pdd_6[k];

        t_12[k] = f_4 * pc_z[k] * pdd_6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_y, pc_x, pc_y, sdd_9, ppf0_5, ppd_2, ppd_9, \
                         ppf1_5, pdd_8, pdd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_0 * sdd_9[k]
                  + f_0 * ppd_9[k]
                  + f_4 * pc_x[k] * pdd_9[k];

        t_14[k] = f_0 * ppd_2[k]
                  + f_4 * pc_y[k] * pdd_8[k];

        t_15[k] = pb_y[k] * ppf0_5[k]
                  - f_5 * pc_y[k] * ppf1_5[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pc_y, pc_z, ppd_3, ppd_5, pdp0_4, pdp0_5, \
                         pdp1_4, pdp1_5, pdd_9, pdd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * ppd_3[k]
                  + f_2 * pdp0_4[k]
                  - f_3 * pdp1_4[k]
                  + f_4 * pc_y[k] * pdd_9[k];

        t_17[k] = f_4 * pc_z[k] * pdd_9[k];

        t_18[k] = f_0 * ppd_5[k]
                  + f_4 * pc_y[k] * pdd_11[k];

        t_19[k] = f_2 * pdp0_5[k]
                  - f_3 * pdp1_5[k]
                  + f_4 * pc_z[k] * pdd_11[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pb_z, pc_y, pc_z, ppf0_0, ppf0_3, \
                         ppd_0, ppf1_0, ppf1_3, pdd_12, pdd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pb_z[k] * ppf0_0[k]
                  - f_5 * pc_z[k] * ppf1_0[k];

        t_21[k] = f_4 * pc_y[k] * pdd_12[k];

        t_22[k] = f_0 * ppd_0[k]
                  + f_4 * pc_z[k] * pdd_12[k];

        t_23[k] = pb_z[k] * ppf0_3[k]
                  - f_5 * pc_z[k] * ppf1_3[k];

        t_24[k] = f_4 * pc_y[k] * pdd_14[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pc_x, pc_y, pc_z, sdd_17, ppd_3, ppd_17, \
                         pdp0_7, pdp1_7, pdd_15, pdd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_0 * sdd_17[k]
                  + f_0 * ppd_17[k]
                  + f_4 * pc_x[k] * pdd_17[k];

        t_26[k] = f_2 * pdp0_7[k]
                  - f_3 * pdp1_7[k]
                  + f_4 * pc_y[k] * pdd_15[k];

        t_27[k] = f_0 * ppd_3[k]
                  + f_4 * pc_z[k] * pdd_15[k];

        t_28[k] = f_4 * pc_y[k] * pdd_17[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pa_x, pc_x, pc_y, pc_z, sdf0_30, sdd_18, sdf1_30, \
                         ppd_5, ppd_6, pdp0_8, pdp1_8, pdd_17, pdd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * ppd_5[k]
                  + f_2 * pdp0_8[k]
                  - f_3 * pdp1_8[k]
                  + f_4 * pc_z[k] * pdd_17[k];

        t_30[k] = pa_x[k] * sdf0_30[k]
                  + f_6 * sdd_18[k]
                  - f_5 * pc_x[k] * sdf1_30[k];

        t_31[k] = f_1 * ppd_6[k]
                  + f_4 * pc_y[k] * pdd_18[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pc_x, pc_y, pc_z, sdd_21, sdd_23, ppd_8, \
                         pdd_18, pdd_20, pdd_21, pdd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_4 * pc_z[k] * pdd_18[k];

        t_33[k] = f_0 * sdd_21[k]
                  + f_4 * pc_x[k] * pdd_21[k];

        t_34[k] = f_1 * ppd_8[k]
                  + f_4 * pc_y[k] * pdd_20[k];

        t_35[k] = f_0 * sdd_23[k]
                  + f_4 * pc_x[k] * pdd_23[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_x, pc_x, pc_y, pc_z, sdf0_36, sdf0_39, \
                         sdf1_36, sdf1_39, ppd_11, pdd_21, pdd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pa_x[k] * sdf0_36[k]
                  - f_5 * pc_x[k] * sdf1_36[k];

        t_37[k] = f_4 * pc_z[k] * pdd_21[k];

        t_38[k] = f_1 * ppd_11[k]
                  + f_4 * pc_y[k] * pdd_23[k];

        t_39[k] = pa_x[k] * sdf0_39[k]
                  - f_5 * pc_x[k] * sdf1_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pb_y, pc_x, pc_y, pc_z, sdd_27, ppf0_20, \
                         ppd_6, ppd_12, ppf1_20, pdd_24, pdd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_y[k] * ppf0_20[k]
                  - f_5 * pc_y[k] * ppf1_20[k];

        t_41[k] = f_0 * ppd_12[k]
                  + f_4 * pc_y[k] * pdd_24[k];

        t_42[k] = f_0 * ppd_6[k]
                  + f_4 * pc_z[k] * pdd_24[k];

        t_43[k] = f_0 * sdd_27[k]
                  + f_4 * pc_x[k] * pdd_27[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pa_x, pc_x, pc_y, pc_z, sdf0_46, sdd_29, \
                         sdf1_46, ppd_9, ppd_14, pdd_26, pdd_27, \
                         pdd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_0 * ppd_14[k]
                  + f_4 * pc_y[k] * pdd_26[k];

        t_45[k] = f_0 * sdd_29[k]
                  + f_4 * pc_x[k] * pdd_29[k];

        t_46[k] = pa_x[k] * sdf0_46[k]
                  - f_5 * pc_x[k] * sdf1_46[k];

        t_47[k] = f_0 * ppd_9[k]
                  + f_4 * pc_z[k] * pdd_27[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_x, pc_x, pc_y, sdf0_49, sdf0_50, sdd_30, \
                         sdf1_49, sdf1_50, ppd_17, pdd_29, pdd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_0 * ppd_17[k]
                  + f_4 * pc_y[k] * pdd_29[k];

        t_49[k] = pa_x[k] * sdf0_49[k]
                  - f_5 * pc_x[k] * sdf1_49[k];

        t_50[k] = pa_x[k] * sdf0_50[k]
                  + f_6 * sdd_30[k]
                  - f_5 * pc_x[k] * sdf1_50[k];

        t_51[k] = f_4 * pc_y[k] * pdd_30[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pc_x, pc_y, pc_z, sdd_33, sdd_35, ppd_12, \
                         pdd_30, pdd_32, pdd_33, pdd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_1 * ppd_12[k]
                  + f_4 * pc_z[k] * pdd_30[k];

        t_53[k] = f_0 * sdd_33[k]
                  + f_4 * pc_x[k] * pdd_33[k];

        t_54[k] = f_4 * pc_y[k] * pdd_32[k];

        t_55[k] = f_0 * sdd_35[k]
                  + f_4 * pc_x[k] * pdd_35[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pa_x, pc_x, pc_y, pc_z, sdf0_56, sdf0_59, \
                         sdf1_56, sdf1_59, ppd_15, pdd_33, pdd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pa_x[k] * sdf0_56[k]
                  - f_5 * pc_x[k] * sdf1_56[k];

        t_57[k] = f_1 * ppd_15[k]
                  + f_4 * pc_z[k] * pdd_33[k];

        t_58[k] = f_4 * pc_y[k] * pdd_35[k];

        t_59[k] = pa_x[k] * sdf0_59[k]
                  - f_5 * pc_x[k] * sdf1_59[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pa_y, pc_x, pc_y, pc_z, sdf0_0, sdf0_1, \
                         sdd_0, sdf1_0, sdf1_1, ppd_21, pdd_36, \
                         pdd_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = pa_y[k] * sdf0_0[k]
                  - f_5 * pc_y[k] * sdf1_0[k];

        t_61[k] = pa_y[k] * sdf0_1[k]
                  + f_0 * sdd_0[k]
                  - f_5 * pc_y[k] * sdf1_1[k];

        t_62[k] = f_4 * pc_z[k] * pdd_36[k];

        t_63[k] = f_1 * ppd_21[k]
                  + f_4 * pc_x[k] * pdd_39[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pa_y, pc_x, pc_y, pc_z, sdf0_6, sdd_3, \
                         sdf1_6, ppd_22, ppd_23, pdd_39, pdd_40, \
                         pdd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_1 * ppd_22[k]
                  + f_4 * pc_x[k] * pdd_40[k];

        t_65[k] = f_1 * ppd_23[k]
                  + f_4 * pc_x[k] * pdd_41[k];

        t_66[k] = pa_y[k] * sdf0_6[k]
                  + f_6 * sdd_3[k]
                  - f_5 * pc_y[k] * sdf1_6[k];

        t_67[k] = f_4 * pc_z[k] * pdd_39[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, pa_y, pc_x, pc_y, sdf0_9, sdd_5, sdf1_9, ppd_24, \
                         pdp0_21, pdp1_21, pdd_41, pdd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_0 * sdd_5[k]
                  + f_4 * pc_y[k] * pdd_41[k];

        t_69[k] = pa_y[k] * sdf0_9[k]
                  - f_5 * pc_y[k] * sdf1_9[k];

        t_70[k] = f_0 * ppd_24[k]
                  + f_2 * pdp0_21[k]
                  - f_3 * pdp1_21[k]
                  + f_4 * pc_x[k] * pdd_42[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, t_74, pb_x, pc_x, pc_z, ppf0_41, ppd_25, ppd_27, \
                         ppd_28, ppf1_41, pdd_42, pdd_45, pdd_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = pb_x[k] * ppf0_41[k]
                  + f_1 * ppd_25[k]
                  - f_5 * pc_x[k] * ppf1_41[k];

        t_72[k] = f_4 * pc_z[k] * pdd_42[k];

        t_73[k] = f_0 * ppd_27[k]
                  + f_4 * pc_x[k] * pdd_45[k];

        t_74[k] = f_0 * ppd_28[k]
                  + f_4 * pc_x[k] * pdd_46[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pb_x, pc_x, pc_z, ppf0_46, ppf0_48, ppd_29, \
                         ppf1_46, ppf1_48, pdd_45, pdd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_0 * ppd_29[k]
                  + f_4 * pc_x[k] * pdd_47[k];

        t_76[k] = pb_x[k] * ppf0_46[k]
                  - f_5 * pc_x[k] * ppf1_46[k];

        t_77[k] = f_4 * pc_z[k] * pdd_45[k];

        t_78[k] = pb_x[k] * ppf0_48[k]
                  - f_5 * pc_x[k] * ppf1_48[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, pa_y, pb_x, pb_z, pc_x, pc_y, pc_z, sdf0_20, \
                         sdf1_20, ppf0_31, ppf0_49, ppf1_31, ppf1_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = pb_x[k] * ppf0_49[k]
                  - f_5 * pc_x[k] * ppf1_49[k];

        t_80[k] = pa_y[k] * sdf0_20[k]
                  - f_5 * pc_y[k] * sdf1_20[k];

        t_81[k] = pb_z[k] * ppf0_31[k]
                  - f_5 * pc_z[k] * ppf1_31[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pc_x, pc_z, ppd_18, ppd_33, ppd_34, ppd_35, \
                         pdd_48, pdd_51, pdd_52, pdd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_0 * ppd_18[k]
                  + f_4 * pc_z[k] * pdd_48[k];

        t_83[k] = f_0 * ppd_33[k]
                  + f_4 * pc_x[k] * pdd_51[k];

        t_84[k] = f_0 * ppd_34[k]
                  + f_4 * pc_x[k] * pdd_52[k];

        t_85[k] = f_0 * ppd_35[k]
                  + f_4 * pc_x[k] * pdd_53[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, pb_x, pc_x, pc_y, pc_z, sdd_17, ppf0_56, ppd_21, \
                         ppf1_56, pdd_51, pdd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pb_x[k] * ppf0_56[k]
                  - f_5 * pc_x[k] * ppf1_56[k];

        t_87[k] = f_0 * ppd_21[k]
                  + f_4 * pc_z[k] * pdd_51[k];

        t_88[k] = f_0 * sdd_17[k]
                  + f_4 * pc_y[k] * pdd_53[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pa_y, pc_x, pc_y, pc_z, sdf0_29, sdf1_29, \
                         pdp0_27, pdp0_28, pdp1_27, pdp1_28, pdd_54, \
                         pdd_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = pa_y[k] * sdf0_29[k]
                  - f_5 * pc_y[k] * sdf1_29[k];

        t_90[k] = f_2 * pdp0_27[k]
                  - f_3 * pdp1_27[k]
                  + f_4 * pc_x[k] * pdd_54[k];

        t_91[k] = f_7 * pdp0_28[k]
                  - f_8 * pdp1_28[k]
                  + f_4 * pc_x[k] * pdd_55[k];

        t_92[k] = f_4 * pc_z[k] * pdd_54[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, pc_x, pc_y, pc_z, sdd_21, ppd_27, \
                         pdp0_28, pdp1_28, pdd_57, pdd_58, pdd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_4 * pc_x[k] * pdd_57[k];

        t_94[k] = f_4 * pc_x[k] * pdd_58[k];

        t_95[k] = f_4 * pc_x[k] * pdd_59[k];

        t_96[k] = f_0 * sdd_21[k]
                  + f_1 * ppd_27[k]
                  + f_2 * pdp0_28[k]
                  - f_3 * pdp1_28[k]
                  + f_4 * pc_y[k] * pdd_57[k];

        t_97[k] = f_4 * pc_z[k] * pdd_57[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, pc_x, pc_y, pc_z, sdd_23, ppd_29, pdp0_29, \
                         pdp0_30, pdp1_29, pdp1_30, pdd_59, pdd_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_0 * sdd_23[k]
                  + f_1 * ppd_29[k]
                  + f_4 * pc_y[k] * pdd_59[k];

        t_99[k] = f_2 * pdp0_29[k]
                  - f_3 * pdp1_29[k]
                  + f_4 * pc_z[k] * pdd_59[k];

        t_100[k] = f_2 * pdp0_30[k]
                   - f_3 * pdp1_30[k]
                   + f_4 * pc_x[k] * pdd_60[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, pb_z, pc_x, pc_z, ppf0_41, ppd_24, \
                         ppf1_41, pdd_60, pdd_63, pdd_64, pdd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = pb_z[k] * ppf0_41[k]
                   - f_5 * pc_z[k] * ppf1_41[k];

        t_102[k] = f_0 * ppd_24[k]
                   + f_4 * pc_z[k] * pdd_60[k];

        t_103[k] = f_4 * pc_x[k] * pdd_63[k];

        t_104[k] = f_4 * pc_x[k] * pdd_64[k];

        t_105[k] = f_4 * pc_x[k] * pdd_65[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, pb_z, pc_y, pc_z, sdd_29, ppf0_46, ppd_27, \
                         ppd_35, ppf1_46, pdd_63, pdd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = pb_z[k] * ppf0_46[k]
                   - f_5 * pc_z[k] * ppf1_46[k];

        t_107[k] = f_0 * ppd_27[k]
                   + f_4 * pc_z[k] * pdd_63[k];

        t_108[k] = f_0 * sdd_29[k]
                   + f_0 * ppd_35[k]
                   + f_4 * pc_y[k] * pdd_65[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, pa_y, pc_x, pc_y, pc_z, sdf0_50, sdf1_50, \
                         ppd_29, pdp0_32, pdp0_34, pdp1_32, pdp1_34, pdd_65, \
                         pdd_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_0 * ppd_29[k]
                   + f_2 * pdp0_32[k]
                   - f_3 * pdp1_32[k]
                   + f_4 * pc_z[k] * pdd_65[k];

        t_110[k] = pa_y[k] * sdf0_50[k]
                   - f_5 * pc_y[k] * sdf1_50[k];

        t_111[k] = f_7 * pdp0_34[k]
                   - f_8 * pdp1_34[k]
                   + f_4 * pc_x[k] * pdd_67[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, t_115, pc_x, pc_z, ppd_30, pdd_66, pdd_69, \
                         pdd_70, pdd_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = f_1 * ppd_30[k]
                   + f_4 * pc_z[k] * pdd_66[k];

        t_113[k] = f_4 * pc_x[k] * pdd_69[k];

        t_114[k] = f_4 * pc_x[k] * pdd_70[k];

        t_115[k] = f_4 * pc_x[k] * pdd_71[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pa_y, pc_y, pc_z, sdf0_56, sdf0_59, \
                         sdd_33, sdd_35, sdf1_56, sdf1_59, ppd_33, pdd_69, \
                         pdd_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = pa_y[k] * sdf0_56[k]
                   + f_6 * sdd_33[k]
                   - f_5 * pc_y[k] * sdf1_56[k];

        t_117[k] = f_1 * ppd_33[k]
                   + f_4 * pc_z[k] * pdd_69[k];

        t_118[k] = f_0 * sdd_35[k]
                   + f_4 * pc_y[k] * pdd_71[k];

        t_119[k] = pa_y[k] * sdf0_59[k]
                   - f_5 * pc_y[k] * sdf1_59[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pa_z, pc_x, pc_y, pc_z, sdf0_0, sdf0_2, \
                         sdd_0, sdf1_0, sdf1_2, ppd_39, pdd_72, \
                         pdd_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = pa_z[k] * sdf0_0[k]
                   - f_5 * pc_z[k] * sdf1_0[k];

        t_121[k] = f_4 * pc_y[k] * pdd_72[k];

        t_122[k] = pa_z[k] * sdf0_2[k]
                   + f_0 * sdd_0[k]
                   - f_5 * pc_z[k] * sdf1_2[k];

        t_123[k] = f_1 * ppd_39[k]
                   + f_4 * pc_x[k] * pdd_75[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, pa_z, pc_x, pc_y, pc_z, sdf0_6, sdf1_6, \
                         ppd_40, ppd_41, pdp0_38, pdp1_38, pdd_76, \
                         pdd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_1 * ppd_40[k]
                   + f_4 * pc_x[k] * pdd_76[k];

        t_125[k] = f_1 * ppd_41[k]
                   + f_4 * pc_x[k] * pdd_77[k];

        t_126[k] = pa_z[k] * sdf0_6[k]
                   - f_5 * pc_z[k] * sdf1_6[k];

        t_127[k] = f_7 * pdp0_38[k]
                   - f_8 * pdp1_38[k]
                   + f_4 * pc_y[k] * pdd_76[k];
    }
}

static auto
compute_prim_pdf_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sdf0, const size_t sdd,
                                                          const size_t sdf1, const size_t ppf0,
                                                          const size_t ppd, const size_t ppf1,
                                                          const size_t pdp0, const size_t pdp1,
                                                          const size_t pdd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.0 / q;
    const auto f_2 = 1.0 / gamma;
    const auto f_3 = p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = gamma / q;
    const auto f_6 = 1.5 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sdf0_9 = buffer.data(sdf0 + 9);
    const auto *sdf0_10 = buffer.data(sdf0 + 10);
    const auto *sdf0_16 = buffer.data(sdf0 + 16);
    const auto *sdf0_30 = buffer.data(sdf0 + 30);
    const auto *sdf0_36 = buffer.data(sdf0 + 36);
    const auto *sdf0_37 = buffer.data(sdf0 + 37);
    const auto *sdf0_39 = buffer.data(sdf0 + 39);

    const auto *sdd_5 = buffer.data(sdd + 5);
    const auto *sdd_21 = buffer.data(sdd + 21);
    const auto *sdd_23 = buffer.data(sdd + 23);
    const auto *sdd_35 = buffer.data(sdd + 35);

    const auto *sdf1_9 = buffer.data(sdf1 + 9);
    const auto *sdf1_10 = buffer.data(sdf1 + 10);
    const auto *sdf1_16 = buffer.data(sdf1 + 16);
    const auto *sdf1_30 = buffer.data(sdf1 + 30);
    const auto *sdf1_36 = buffer.data(sdf1 + 36);
    const auto *sdf1_37 = buffer.data(sdf1 + 37);
    const auto *sdf1_39 = buffer.data(sdf1 + 39);

    const auto *ppf0_62 = buffer.data(ppf0 + 62);
    const auto *ppf0_77 = buffer.data(ppf0 + 77);
    const auto *ppf0_79 = buffer.data(ppf0 + 79);
    const auto *ppf0_80 = buffer.data(ppf0 + 80);
    const auto *ppf0_82 = buffer.data(ppf0 + 82);
    const auto *ppf0_86 = buffer.data(ppf0 + 86);
    const auto *ppf0_87 = buffer.data(ppf0 + 87);
    const auto *ppf0_89 = buffer.data(ppf0 + 89);

    const auto *ppd_36 = buffer.data(ppd + 36);
    const auto *ppd_41 = buffer.data(ppd + 41);
    const auto *ppd_42 = buffer.data(ppd + 42);
    const auto *ppd_45 = buffer.data(ppd + 45);
    const auto *ppd_46 = buffer.data(ppd + 46);
    const auto *ppd_47 = buffer.data(ppd + 47);
    const auto *ppd_48 = buffer.data(ppd + 48);
    const auto *ppd_50 = buffer.data(ppd + 50);
    const auto *ppd_51 = buffer.data(ppd + 51);
    const auto *ppd_52 = buffer.data(ppd + 52);
    const auto *ppd_53 = buffer.data(ppd + 53);

    const auto *ppf1_62 = buffer.data(ppf1 + 62);
    const auto *ppf1_77 = buffer.data(ppf1 + 77);
    const auto *ppf1_79 = buffer.data(ppf1 + 79);
    const auto *ppf1_80 = buffer.data(ppf1 + 80);
    const auto *ppf1_82 = buffer.data(ppf1 + 82);
    const auto *ppf1_86 = buffer.data(ppf1 + 86);
    const auto *ppf1_87 = buffer.data(ppf1 + 87);
    const auto *ppf1_89 = buffer.data(ppf1 + 89);

    const auto *pdp0_42 = buffer.data(pdp0 + 42);
    const auto *pdp0_47 = buffer.data(pdp0 + 47);
    const auto *pdp0_49 = buffer.data(pdp0 + 49);
    const auto *pdp0_51 = buffer.data(pdp0 + 51);
    const auto *pdp0_52 = buffer.data(pdp0 + 52);
    const auto *pdp0_53 = buffer.data(pdp0 + 53);

    const auto *pdp1_42 = buffer.data(pdp1 + 42);
    const auto *pdp1_47 = buffer.data(pdp1 + 47);
    const auto *pdp1_49 = buffer.data(pdp1 + 49);
    const auto *pdp1_51 = buffer.data(pdp1 + 51);
    const auto *pdp1_52 = buffer.data(pdp1 + 52);
    const auto *pdp1_53 = buffer.data(pdp1 + 53);

    const auto *pdd_77 = buffer.data(pdd + 77);
    const auto *pdd_78 = buffer.data(pdd + 78);
    const auto *pdd_81 = buffer.data(pdd + 81);
    const auto *pdd_82 = buffer.data(pdd + 82);
    const auto *pdd_83 = buffer.data(pdd + 83);
    const auto *pdd_84 = buffer.data(pdd + 84);
    const auto *pdd_87 = buffer.data(pdd + 87);
    const auto *pdd_88 = buffer.data(pdd + 88);
    const auto *pdd_89 = buffer.data(pdd + 89);
    const auto *pdd_90 = buffer.data(pdd + 90);
    const auto *pdd_92 = buffer.data(pdd + 92);
    const auto *pdd_93 = buffer.data(pdd + 93);
    const auto *pdd_94 = buffer.data(pdd + 94);
    const auto *pdd_95 = buffer.data(pdd + 95);
    const auto *pdd_96 = buffer.data(pdd + 96);
    const auto *pdd_99 = buffer.data(pdd + 99);
    const auto *pdd_100 = buffer.data(pdd + 100);
    const auto *pdd_101 = buffer.data(pdd + 101);
    const auto *pdd_102 = buffer.data(pdd + 102);
    const auto *pdd_104 = buffer.data(pdd + 104);
    const auto *pdd_105 = buffer.data(pdd + 105);
    const auto *pdd_106 = buffer.data(pdd + 106);
    const auto *pdd_107 = buffer.data(pdd + 107);

#pragma omp simd aligned(t_128, t_129, t_130, t_131, pa_z, pc_y, pc_z, sdf0_9, sdf0_10, sdd_5, \
                         sdf1_9, sdf1_10, ppd_36, pdd_77, pdd_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_4 * pc_y[k] * pdd_77[k];

        t_129[k] = pa_z[k] * sdf0_9[k]
                   + f_6 * sdd_5[k]
                   - f_5 * pc_z[k] * sdf1_9[k];

        t_130[k] = pa_z[k] * sdf0_10[k]
                   - f_5 * pc_z[k] * sdf1_10[k];

        t_131[k] = f_0 * ppd_36[k]
                   + f_4 * pc_y[k] * pdd_78[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pb_y, pc_x, pc_y, ppf0_62, ppd_45, \
                         ppd_46, ppd_47, ppf1_62, pdd_81, pdd_82, \
                         pdd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = pb_y[k] * ppf0_62[k]
                   - f_5 * pc_y[k] * ppf1_62[k];

        t_133[k] = f_0 * ppd_45[k]
                   + f_4 * pc_x[k] * pdd_81[k];

        t_134[k] = f_0 * ppd_46[k]
                   + f_4 * pc_x[k] * pdd_82[k];

        t_135[k] = f_0 * ppd_47[k]
                   + f_4 * pc_x[k] * pdd_83[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pa_z, pb_x, pc_x, pc_y, pc_z, sdf0_16, sdf1_16, \
                         ppf0_77, ppd_41, ppf1_77, pdd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = pa_z[k] * sdf0_16[k]
                   - f_5 * pc_z[k] * sdf1_16[k];

        t_137[k] = pb_x[k] * ppf0_77[k]
                   - f_5 * pc_x[k] * ppf1_77[k];

        t_138[k] = f_0 * ppd_41[k]
                   + f_4 * pc_y[k] * pdd_83[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pb_x, pc_x, pc_y, ppf0_79, ppf0_82, \
                         ppd_48, ppd_50, ppf1_79, ppf1_82, pdp0_42, pdp1_42, \
                         pdd_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pb_x[k] * ppf0_79[k]
                   - f_5 * pc_x[k] * ppf1_79[k];

        t_140[k] = f_0 * ppd_48[k]
                   + f_2 * pdp0_42[k]
                   - f_3 * pdp1_42[k]
                   + f_4 * pc_x[k] * pdd_84[k];

        t_141[k] = f_4 * pc_y[k] * pdd_84[k];

        t_142[k] = pb_x[k] * ppf0_82[k]
                   + f_1 * ppd_50[k]
                   - f_5 * pc_x[k] * ppf1_82[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pb_x, pc_x, ppf0_86, ppd_51, ppd_52, \
                         ppd_53, ppf1_86, pdd_87, pdd_88, pdd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_0 * ppd_51[k]
                   + f_4 * pc_x[k] * pdd_87[k];

        t_144[k] = f_0 * ppd_52[k]
                   + f_4 * pc_x[k] * pdd_88[k];

        t_145[k] = f_0 * ppd_53[k]
                   + f_4 * pc_x[k] * pdd_89[k];

        t_146[k] = pb_x[k] * ppf0_86[k]
                   - f_5 * pc_x[k] * ppf1_86[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, pa_z, pb_x, pc_x, pc_y, pc_z, sdf0_30, \
                         sdf1_30, ppf0_87, ppf0_89, ppf1_87, ppf1_89, \
                         pdd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = pb_x[k] * ppf0_87[k]
                   - f_5 * pc_x[k] * ppf1_87[k];

        t_148[k] = f_4 * pc_y[k] * pdd_89[k];

        t_149[k] = pb_x[k] * ppf0_89[k]
                   - f_5 * pc_x[k] * ppf1_89[k];

        t_150[k] = pa_z[k] * sdf0_30[k]
                   - f_5 * pc_z[k] * sdf1_30[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, t_155, pc_x, pc_y, ppd_42, pdp0_47, \
                         pdp1_47, pdd_90, pdd_92, pdd_93, pdd_94, \
                         pdd_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_1 * ppd_42[k]
                   + f_4 * pc_y[k] * pdd_90[k];

        t_152[k] = f_7 * pdp0_47[k]
                   - f_8 * pdp1_47[k]
                   + f_4 * pc_x[k] * pdd_92[k];

        t_153[k] = f_4 * pc_x[k] * pdd_93[k];

        t_154[k] = f_4 * pc_x[k] * pdd_94[k];

        t_155[k] = f_4 * pc_x[k] * pdd_95[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pa_z, pc_y, pc_z, sdf0_36, sdf0_37, sdd_21, \
                         sdf1_36, sdf1_37, ppd_47, pdd_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = pa_z[k] * sdf0_36[k]
                   - f_5 * pc_z[k] * sdf1_36[k];

        t_157[k] = pa_z[k] * sdf0_37[k]
                   + f_0 * sdd_21[k]
                   - f_5 * pc_z[k] * sdf1_37[k];

        t_158[k] = f_1 * ppd_47[k]
                   + f_4 * pc_y[k] * pdd_95[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, pa_z, pb_y, pc_y, pc_z, sdf0_39, sdd_23, \
                         sdf1_39, ppf0_80, ppd_48, ppf1_80, pdd_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = pa_z[k] * sdf0_39[k]
                   + f_6 * sdd_23[k]
                   - f_5 * pc_z[k] * sdf1_39[k];

        t_160[k] = pb_y[k] * ppf0_80[k]
                   - f_5 * pc_y[k] * ppf1_80[k];

        t_161[k] = f_0 * ppd_48[k]
                   + f_4 * pc_y[k] * pdd_96[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, pb_y, pc_x, pc_y, ppf0_82, ppd_51, \
                         ppf1_82, pdp0_49, pdp1_49, pdd_99, pdd_100, \
                         pdd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = pb_y[k] * ppf0_82[k]
                   - f_5 * pc_y[k] * ppf1_82[k];

        t_163[k] = f_4 * pc_x[k] * pdd_99[k];

        t_164[k] = f_4 * pc_x[k] * pdd_100[k];

        t_165[k] = f_4 * pc_x[k] * pdd_101[k];

        t_166[k] = f_0 * ppd_51[k]
                   + f_2 * pdp0_49[k]
                   - f_3 * pdp1_49[k]
                   + f_4 * pc_y[k] * pdd_99[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, pb_y, pc_y, ppf0_87, ppf0_89, ppd_52, ppd_53, \
                         ppf1_87, ppf1_89, pdd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = pb_y[k] * ppf0_87[k]
                   + f_1 * ppd_52[k]
                   - f_5 * pc_y[k] * ppf1_87[k];

        t_168[k] = f_0 * ppd_53[k]
                   + f_4 * pc_y[k] * pdd_101[k];

        t_169[k] = pb_y[k] * ppf0_89[k]
                   - f_5 * pc_y[k] * ppf1_89[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, pc_x, pc_y, pdp0_51, pdp0_53, \
                         pdp1_51, pdp1_53, pdd_102, pdd_104, pdd_105, \
                         pdd_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_2 * pdp0_51[k]
                   - f_3 * pdp1_51[k]
                   + f_4 * pc_x[k] * pdd_102[k];

        t_171[k] = f_4 * pc_y[k] * pdd_102[k];

        t_172[k] = f_7 * pdp0_53[k]
                   - f_8 * pdp1_53[k]
                   + f_4 * pc_x[k] * pdd_104[k];

        t_173[k] = f_4 * pc_x[k] * pdd_105[k];

        t_174[k] = f_4 * pc_x[k] * pdd_106[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, pc_x, pc_y, pdp0_52, pdp0_53, pdp1_52, \
                         pdp1_53, pdd_105, pdd_106, pdd_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_4 * pc_x[k] * pdd_107[k];

        t_176[k] = f_2 * pdp0_52[k]
                   - f_3 * pdp1_52[k]
                   + f_4 * pc_y[k] * pdd_105[k];

        t_177[k] = f_7 * pdp0_53[k]
                   - f_8 * pdp1_53[k]
                   + f_4 * pc_y[k] * pdd_106[k];

        t_178[k] = f_4 * pc_y[k] * pdd_107[k];
    }

#pragma omp simd aligned(t_179, pc_z, sdd_35, ppd_53, pdp0_53, pdp1_53, \
                         pdd_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_0 * sdd_35[k]
                   + f_1 * ppd_53[k]
                   + f_2 * pdp0_53[k]
                   - f_3 * pdp1_53[k]
                   + f_4 * pc_z[k] * pdd_107[k];
    }
}

auto
compute_prim_pdf_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t sdf0,
                                                   const size_t sdd, const size_t sdf1,
                                                   const size_t ppf0, const size_t ppd,
                                                   const size_t ppf1, const size_t pdp0,
                                                   const size_t pdp1, const size_t pdd,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_pdf_three_center_electron_repulsion_0_piece0(buffer, target, pa, pb, pc, sdf0,
                                                              sdd, sdf1, ppf0, ppd, ppf1, pdp0,
                                                              pdp1, pdd, ncols, gamma, p, q);

    compute_prim_pdf_three_center_electron_repulsion_0_piece1(buffer, target, pa, pb, pc, sdf0,
                                                              sdd, sdf1, ppf0, ppd, ppf1, pdp0,
                                                              pdp1, pdd, ncols, gamma, p, q);
}

}  // namespace simdt3ceri
