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


#include "SimdThreeCenterElectronRepulsionVrrRecPDD.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_pdd_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t sdd0,
                                                   const size_t sdp, const size_t sdd1,
                                                   const size_t ppd0, const size_t ppp,
                                                   const size_t ppd1, const size_t pds0,
                                                   const size_t pds1, const size_t pdp,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.0 / q;
    const auto f_2 = 0.5 / gamma;
    const auto f_3 = 0.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = gamma / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sdd0_0 = buffer.data(sdd0 + 0);
    const auto *sdd0_3 = buffer.data(sdd0 + 3);
    const auto *sdd0_5 = buffer.data(sdd0 + 5);
    const auto *sdd0_6 = buffer.data(sdd0 + 6);
    const auto *sdd0_9 = buffer.data(sdd0 + 9);
    const auto *sdd0_12 = buffer.data(sdd0 + 12);
    const auto *sdd0_17 = buffer.data(sdd0 + 17);
    const auto *sdd0_18 = buffer.data(sdd0 + 18);
    const auto *sdd0_21 = buffer.data(sdd0 + 21);
    const auto *sdd0_23 = buffer.data(sdd0 + 23);
    const auto *sdd0_27 = buffer.data(sdd0 + 27);
    const auto *sdd0_29 = buffer.data(sdd0 + 29);
    const auto *sdd0_30 = buffer.data(sdd0 + 30);
    const auto *sdd0_33 = buffer.data(sdd0 + 33);
    const auto *sdd0_35 = buffer.data(sdd0 + 35);

    const auto *sdp_0 = buffer.data(sdp + 0);
    const auto *sdp_1 = buffer.data(sdp + 1);
    const auto *sdp_2 = buffer.data(sdp + 2);
    const auto *sdp_9 = buffer.data(sdp + 9);
    const auto *sdp_10 = buffer.data(sdp + 10);
    const auto *sdp_11 = buffer.data(sdp + 11);
    const auto *sdp_15 = buffer.data(sdp + 15);
    const auto *sdp_16 = buffer.data(sdp + 16);
    const auto *sdp_17 = buffer.data(sdp + 17);

    const auto *sdd1_0 = buffer.data(sdd1 + 0);
    const auto *sdd1_3 = buffer.data(sdd1 + 3);
    const auto *sdd1_5 = buffer.data(sdd1 + 5);
    const auto *sdd1_6 = buffer.data(sdd1 + 6);
    const auto *sdd1_9 = buffer.data(sdd1 + 9);
    const auto *sdd1_12 = buffer.data(sdd1 + 12);
    const auto *sdd1_17 = buffer.data(sdd1 + 17);
    const auto *sdd1_18 = buffer.data(sdd1 + 18);
    const auto *sdd1_21 = buffer.data(sdd1 + 21);
    const auto *sdd1_23 = buffer.data(sdd1 + 23);
    const auto *sdd1_27 = buffer.data(sdd1 + 27);
    const auto *sdd1_29 = buffer.data(sdd1 + 29);
    const auto *sdd1_30 = buffer.data(sdd1 + 30);
    const auto *sdd1_33 = buffer.data(sdd1 + 33);
    const auto *sdd1_35 = buffer.data(sdd1 + 35);

    const auto *ppd0_0 = buffer.data(ppd0 + 0);
    const auto *ppd0_12 = buffer.data(ppd0 + 12);
    const auto *ppd0_27 = buffer.data(ppd0 + 27);
    const auto *ppd0_29 = buffer.data(ppd0 + 29);
    const auto *ppd0_33 = buffer.data(ppd0 + 33);
    const auto *ppd0_47 = buffer.data(ppd0 + 47);
    const auto *ppd0_48 = buffer.data(ppd0 + 48);
    const auto *ppd0_51 = buffer.data(ppd0 + 51);
    const auto *ppd0_53 = buffer.data(ppd0 + 53);

    const auto *ppp_0 = buffer.data(ppp + 0);
    const auto *ppp_1 = buffer.data(ppp + 1);
    const auto *ppp_2 = buffer.data(ppp + 2);
    const auto *ppp_3 = buffer.data(ppp + 3);
    const auto *ppp_5 = buffer.data(ppp + 5);
    const auto *ppp_6 = buffer.data(ppp + 6);
    const auto *ppp_8 = buffer.data(ppp + 8);
    const auto *ppp_10 = buffer.data(ppp + 10);
    const auto *ppp_11 = buffer.data(ppp + 11);
    const auto *ppp_12 = buffer.data(ppp + 12);
    const auto *ppp_13 = buffer.data(ppp + 13);
    const auto *ppp_14 = buffer.data(ppp + 14);
    const auto *ppp_16 = buffer.data(ppp + 16);
    const auto *ppp_17 = buffer.data(ppp + 17);
    const auto *ppp_19 = buffer.data(ppp + 19);
    const auto *ppp_20 = buffer.data(ppp + 20);
    const auto *ppp_22 = buffer.data(ppp + 22);
    const auto *ppp_23 = buffer.data(ppp + 23);
    const auto *ppp_24 = buffer.data(ppp + 24);
    const auto *ppp_25 = buffer.data(ppp + 25);
    const auto *ppp_26 = buffer.data(ppp + 26);

    const auto *ppd1_0 = buffer.data(ppd1 + 0);
    const auto *ppd1_12 = buffer.data(ppd1 + 12);
    const auto *ppd1_27 = buffer.data(ppd1 + 27);
    const auto *ppd1_29 = buffer.data(ppd1 + 29);
    const auto *ppd1_33 = buffer.data(ppd1 + 33);
    const auto *ppd1_47 = buffer.data(ppd1 + 47);
    const auto *ppd1_48 = buffer.data(ppd1 + 48);
    const auto *ppd1_51 = buffer.data(ppd1 + 51);
    const auto *ppd1_53 = buffer.data(ppd1 + 53);

    const auto *pds0_0 = buffer.data(pds0 + 0);
    const auto *pds0_1 = buffer.data(pds0 + 1);
    const auto *pds0_2 = buffer.data(pds0 + 2);
    const auto *pds0_7 = buffer.data(pds0 + 7);
    const auto *pds0_9 = buffer.data(pds0 + 9);
    const auto *pds0_10 = buffer.data(pds0 + 10);
    const auto *pds0_14 = buffer.data(pds0 + 14);
    const auto *pds0_16 = buffer.data(pds0 + 16);
    const auto *pds0_17 = buffer.data(pds0 + 17);

    const auto *pds1_0 = buffer.data(pds1 + 0);
    const auto *pds1_1 = buffer.data(pds1 + 1);
    const auto *pds1_2 = buffer.data(pds1 + 2);
    const auto *pds1_7 = buffer.data(pds1 + 7);
    const auto *pds1_9 = buffer.data(pds1 + 9);
    const auto *pds1_10 = buffer.data(pds1 + 10);
    const auto *pds1_14 = buffer.data(pds1 + 14);
    const auto *pds1_16 = buffer.data(pds1 + 16);
    const auto *pds1_17 = buffer.data(pds1 + 17);

    const auto *pdp_0 = buffer.data(pdp + 0);
    const auto *pdp_1 = buffer.data(pdp + 1);
    const auto *pdp_2 = buffer.data(pdp + 2);
    const auto *pdp_3 = buffer.data(pdp + 3);
    const auto *pdp_4 = buffer.data(pdp + 4);
    const auto *pdp_5 = buffer.data(pdp + 5);
    const auto *pdp_6 = buffer.data(pdp + 6);
    const auto *pdp_7 = buffer.data(pdp + 7);
    const auto *pdp_8 = buffer.data(pdp + 8);
    const auto *pdp_9 = buffer.data(pdp + 9);
    const auto *pdp_11 = buffer.data(pdp + 11);
    const auto *pdp_12 = buffer.data(pdp + 12);
    const auto *pdp_14 = buffer.data(pdp + 14);
    const auto *pdp_15 = buffer.data(pdp + 15);
    const auto *pdp_17 = buffer.data(pdp + 17);
    const auto *pdp_19 = buffer.data(pdp + 19);
    const auto *pdp_20 = buffer.data(pdp + 20);
    const auto *pdp_21 = buffer.data(pdp + 21);
    const auto *pdp_22 = buffer.data(pdp + 22);
    const auto *pdp_23 = buffer.data(pdp + 23);
    const auto *pdp_25 = buffer.data(pdp + 25);
    const auto *pdp_26 = buffer.data(pdp + 26);
    const auto *pdp_27 = buffer.data(pdp + 27);
    const auto *pdp_28 = buffer.data(pdp + 28);
    const auto *pdp_29 = buffer.data(pdp + 29);
    const auto *pdp_30 = buffer.data(pdp + 30);
    const auto *pdp_31 = buffer.data(pdp + 31);
    const auto *pdp_32 = buffer.data(pdp + 32);
    const auto *pdp_34 = buffer.data(pdp + 34);
    const auto *pdp_35 = buffer.data(pdp + 35);
    const auto *pdp_37 = buffer.data(pdp + 37);
    const auto *pdp_38 = buffer.data(pdp + 38);
    const auto *pdp_40 = buffer.data(pdp + 40);
    const auto *pdp_41 = buffer.data(pdp + 41);
    const auto *pdp_42 = buffer.data(pdp + 42);
    const auto *pdp_43 = buffer.data(pdp + 43);
    const auto *pdp_44 = buffer.data(pdp + 44);
    const auto *pdp_46 = buffer.data(pdp + 46);
    const auto *pdp_47 = buffer.data(pdp + 47);
    const auto *pdp_49 = buffer.data(pdp + 49);
    const auto *pdp_50 = buffer.data(pdp + 50);
    const auto *pdp_51 = buffer.data(pdp + 51);
    const auto *pdp_52 = buffer.data(pdp + 52);
    const auto *pdp_53 = buffer.data(pdp + 53);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, sdp_0, ppp_0, pds0_0, \
                         pds1_0, pdp_0, pdp_1, pdp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sdp_0[k]
                 + f_1 * ppp_0[k]
                 + f_2 * pds0_0[k]
                 - f_3 * pds1_0[k]
                 + f_4 * pc_x[k] * pdp_0[k];

        t_1[k] = f_4 * pc_y[k] * pdp_0[k];

        t_2[k] = f_4 * pc_z[k] * pdp_0[k];

        t_3[k] = f_2 * pds0_0[k]
                 - f_3 * pds1_0[k]
                 + f_4 * pc_y[k] * pdp_1[k];

        t_4[k] = f_4 * pc_y[k] * pdp_2[k];

        t_5[k] = f_2 * pds0_0[k]
                 - f_3 * pds1_0[k]
                 + f_4 * pc_z[k] * pdp_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pb_y, pc_y, pc_z, ppd0_0, ppp_0, ppp_1, ppd1_0, \
                         pds0_1, pds1_1, pdp_3, pdp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pb_y[k] * ppd0_0[k]
                 - f_5 * pc_y[k] * ppd1_0[k];

        t_7[k] = f_0 * ppp_0[k]
                 + f_4 * pc_y[k] * pdp_3[k];

        t_8[k] = f_4 * pc_z[k] * pdp_3[k];

        t_9[k] = f_0 * ppp_1[k]
                 + f_2 * pds0_1[k]
                 - f_3 * pds1_1[k]
                 + f_4 * pc_y[k] * pdp_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pb_z, pc_y, pc_z, ppd0_0, ppp_0, ppp_2, \
                         ppd1_0, pds0_1, pds1_1, pdp_5, pdp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * ppp_2[k]
                  + f_4 * pc_y[k] * pdp_5[k];

        t_11[k] = f_2 * pds0_1[k]
                  - f_3 * pds1_1[k]
                  + f_4 * pc_z[k] * pdp_5[k];

        t_12[k] = pb_z[k] * ppd0_0[k]
                  - f_5 * pc_z[k] * ppd1_0[k];

        t_13[k] = f_4 * pc_y[k] * pdp_6[k];

        t_14[k] = f_0 * ppp_0[k]
                  + f_4 * pc_z[k] * pdp_6[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_x, pc_x, pc_y, pc_z, sdd0_18, sdp_9, \
                         sdd1_18, ppp_2, pds0_2, pds1_2, pdp_7, pdp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_2 * pds0_2[k]
                  - f_3 * pds1_2[k]
                  + f_4 * pc_y[k] * pdp_7[k];

        t_16[k] = f_4 * pc_y[k] * pdp_8[k];

        t_17[k] = f_0 * ppp_2[k]
                  + f_2 * pds0_2[k]
                  - f_3 * pds1_2[k]
                  + f_4 * pc_z[k] * pdp_8[k];

        t_18[k] = pa_x[k] * sdd0_18[k]
                  + f_1 * sdp_9[k]
                  - f_5 * pc_x[k] * sdd1_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_x, pc_x, pc_y, pc_z, sdd0_21, sdd1_21, \
                         ppp_3, ppp_5, pdp_9, pdp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_1 * ppp_3[k]
                  + f_4 * pc_y[k] * pdp_9[k];

        t_20[k] = f_4 * pc_z[k] * pdp_9[k];

        t_21[k] = pa_x[k] * sdd0_21[k]
                  - f_5 * pc_x[k] * sdd1_21[k];

        t_22[k] = f_1 * ppp_5[k]
                  + f_4 * pc_y[k] * pdp_11[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pb_y, pc_x, pc_y, pc_z, sdd0_23, \
                         sdd1_23, ppd0_12, ppp_3, ppp_6, ppd1_12, \
                         pdp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = pa_x[k] * sdd0_23[k]
                  - f_5 * pc_x[k] * sdd1_23[k];

        t_24[k] = pb_y[k] * ppd0_12[k]
                  - f_5 * pc_y[k] * ppd1_12[k];

        t_25[k] = f_0 * ppp_6[k]
                  + f_4 * pc_y[k] * pdp_12[k];

        t_26[k] = f_0 * ppp_3[k]
                  + f_4 * pc_z[k] * pdp_12[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pa_x, pc_x, pc_y, sdd0_27, sdd0_29, sdd0_30, \
                         sdp_15, sdd1_27, sdd1_29, sdd1_30, ppp_8, \
                         pdp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pa_x[k] * sdd0_27[k]
                  - f_5 * pc_x[k] * sdd1_27[k];

        t_28[k] = f_0 * ppp_8[k]
                  + f_4 * pc_y[k] * pdp_14[k];

        t_29[k] = pa_x[k] * sdd0_29[k]
                  - f_5 * pc_x[k] * sdd1_29[k];

        t_30[k] = pa_x[k] * sdd0_30[k]
                  + f_1 * sdp_15[k]
                  - f_5 * pc_x[k] * sdd1_30[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pa_x, pc_x, pc_y, pc_z, sdd0_33, \
                         sdd0_35, sdd1_33, sdd1_35, ppp_6, pdp_15, \
                         pdp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_4 * pc_y[k] * pdp_15[k];

        t_32[k] = f_1 * ppp_6[k]
                  + f_4 * pc_z[k] * pdp_15[k];

        t_33[k] = pa_x[k] * sdd0_33[k]
                  - f_5 * pc_x[k] * sdd1_33[k];

        t_34[k] = f_4 * pc_y[k] * pdp_17[k];

        t_35[k] = pa_x[k] * sdd0_35[k]
                  - f_5 * pc_x[k] * sdd1_35[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_y, pc_x, pc_y, sdd0_0, sdd0_3, sdp_1, \
                         sdd1_0, sdd1_3, ppp_10, ppp_11, pdp_19, \
                         pdp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pa_y[k] * sdd0_0[k]
                  - f_5 * pc_y[k] * sdd1_0[k];

        t_37[k] = f_1 * ppp_10[k]
                  + f_4 * pc_x[k] * pdp_19[k];

        t_38[k] = f_1 * ppp_11[k]
                  + f_4 * pc_x[k] * pdp_20[k];

        t_39[k] = pa_y[k] * sdd0_3[k]
                  + f_1 * sdp_1[k]
                  - f_5 * pc_y[k] * sdd1_3[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_y, pc_x, pc_y, pc_z, sdd0_5, sdd1_5, ppp_12, \
                         pds0_7, pds1_7, pdp_19, pdp_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_4 * pc_z[k] * pdp_19[k];

        t_41[k] = pa_y[k] * sdd0_5[k]
                  - f_5 * pc_y[k] * sdd1_5[k];

        t_42[k] = f_0 * ppp_12[k]
                  + f_2 * pds0_7[k]
                  - f_3 * pds1_7[k]
                  + f_4 * pc_x[k] * pdp_21[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, t_47, pb_x, pc_x, pc_z, ppd0_27, ppd0_29, \
                         ppp_13, ppp_14, ppd1_27, ppd1_29, pdp_22, \
                         pdp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_0 * ppp_13[k]
                  + f_4 * pc_x[k] * pdp_22[k];

        t_44[k] = f_0 * ppp_14[k]
                  + f_4 * pc_x[k] * pdp_23[k];

        t_45[k] = pb_x[k] * ppd0_27[k]
                  - f_5 * pc_x[k] * ppd1_27[k];

        t_46[k] = f_4 * pc_z[k] * pdp_22[k];

        t_47[k] = pb_x[k] * ppd0_29[k]
                  - f_5 * pc_x[k] * ppd1_29[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, t_51, pa_y, pb_x, pc_x, pc_y, sdd0_12, sdd1_12, \
                         ppd0_33, ppp_16, ppp_17, ppd1_33, pdp_25, \
                         pdp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = pa_y[k] * sdd0_12[k]
                  - f_5 * pc_y[k] * sdd1_12[k];

        t_49[k] = f_0 * ppp_16[k]
                  + f_4 * pc_x[k] * pdp_25[k];

        t_50[k] = f_0 * ppp_17[k]
                  + f_4 * pc_x[k] * pdp_26[k];

        t_51[k] = pb_x[k] * ppd0_33[k]
                  - f_5 * pc_x[k] * ppd1_33[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pa_y, pc_x, pc_y, pc_z, sdd0_17, sdd1_17, \
                         ppp_10, pds0_9, pds1_9, pdp_25, pdp_27, \
                         pdp_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_0 * ppp_10[k]
                  + f_4 * pc_z[k] * pdp_25[k];

        t_53[k] = pa_y[k] * sdd0_17[k]
                  - f_5 * pc_y[k] * sdd1_17[k];

        t_54[k] = f_2 * pds0_9[k]
                  - f_3 * pds1_9[k]
                  + f_4 * pc_x[k] * pdp_27[k];

        t_55[k] = f_4 * pc_x[k] * pdp_28[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pc_x, pc_y, pc_z, sdp_10, ppp_13, pds0_9, \
                         pds1_9, pdp_28, pdp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_4 * pc_x[k] * pdp_29[k];

        t_57[k] = f_0 * sdp_10[k]
                  + f_1 * ppp_13[k]
                  + f_2 * pds0_9[k]
                  - f_3 * pds1_9[k]
                  + f_4 * pc_y[k] * pdp_28[k];

        t_58[k] = f_4 * pc_z[k] * pdp_28[k];

        t_59[k] = f_2 * pds0_9[k]
                  - f_3 * pds1_9[k]
                  + f_4 * pc_z[k] * pdp_29[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, t_64, pb_z, pc_x, pc_z, ppd0_27, ppp_13, \
                         ppd1_27, pds0_10, pds1_10, pdp_30, pdp_31, \
                         pdp_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_2 * pds0_10[k]
                  - f_3 * pds1_10[k]
                  + f_4 * pc_x[k] * pdp_30[k];

        t_61[k] = f_4 * pc_x[k] * pdp_31[k];

        t_62[k] = f_4 * pc_x[k] * pdp_32[k];

        t_63[k] = pb_z[k] * ppd0_27[k]
                  - f_5 * pc_z[k] * ppd1_27[k];

        t_64[k] = f_0 * ppp_13[k]
                  + f_4 * pc_z[k] * pdp_31[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pa_y, pc_x, pc_y, pc_z, sdd0_30, sdd1_30, \
                         ppp_14, pds0_10, pds1_10, pdp_32, pdp_34, \
                         pdp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_0 * ppp_14[k]
                  + f_2 * pds0_10[k]
                  - f_3 * pds1_10[k]
                  + f_4 * pc_z[k] * pdp_32[k];

        t_66[k] = pa_y[k] * sdd0_30[k]
                  - f_5 * pc_y[k] * sdd1_30[k];

        t_67[k] = f_4 * pc_x[k] * pdp_34[k];

        t_68[k] = f_4 * pc_x[k] * pdp_35[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pa_y, pc_y, pc_z, sdd0_33, sdd0_35, sdp_16, \
                         sdd1_33, sdd1_35, ppp_16, pdp_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = pa_y[k] * sdd0_33[k]
                  + f_1 * sdp_16[k]
                  - f_5 * pc_y[k] * sdd1_33[k];

        t_70[k] = f_1 * ppp_16[k]
                  + f_4 * pc_z[k] * pdp_34[k];

        t_71[k] = pa_y[k] * sdd0_35[k]
                  - f_5 * pc_y[k] * sdd1_35[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pa_z, pc_x, pc_z, sdd0_0, sdd0_3, sdd1_0, \
                         sdd1_3, ppp_19, ppp_20, pdp_37, pdp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = pa_z[k] * sdd0_0[k]
                  - f_5 * pc_z[k] * sdd1_0[k];

        t_73[k] = f_1 * ppp_19[k]
                  + f_4 * pc_x[k] * pdp_37[k];

        t_74[k] = f_1 * ppp_20[k]
                  + f_4 * pc_x[k] * pdp_38[k];

        t_75[k] = pa_z[k] * sdd0_3[k]
                  - f_5 * pc_z[k] * sdd1_3[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pa_z, pc_x, pc_y, pc_z, sdd0_5, sdd0_6, \
                         sdp_2, sdd1_5, sdd1_6, ppp_22, pdp_38, \
                         pdp_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_4 * pc_y[k] * pdp_38[k];

        t_77[k] = pa_z[k] * sdd0_5[k]
                  + f_1 * sdp_2[k]
                  - f_5 * pc_z[k] * sdd1_5[k];

        t_78[k] = pa_z[k] * sdd0_6[k]
                  - f_5 * pc_z[k] * sdd1_6[k];

        t_79[k] = f_0 * ppp_22[k]
                  + f_4 * pc_x[k] * pdp_40[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, pa_z, pb_x, pc_x, pc_y, pc_z, sdd0_9, sdd1_9, \
                         ppd0_47, ppp_20, ppp_23, ppd1_47, pdp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_0 * ppp_23[k]
                  + f_4 * pc_x[k] * pdp_41[k];

        t_81[k] = pa_z[k] * sdd0_9[k]
                  - f_5 * pc_z[k] * sdd1_9[k];

        t_82[k] = f_0 * ppp_20[k]
                  + f_4 * pc_y[k] * pdp_41[k];

        t_83[k] = pb_x[k] * ppd0_47[k]
                  - f_5 * pc_x[k] * ppd1_47[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, pb_x, pc_x, ppd0_51, ppp_24, ppp_25, ppp_26, \
                         ppd1_51, pds0_14, pds1_14, pdp_42, pdp_43, \
                         pdp_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_0 * ppp_24[k]
                  + f_2 * pds0_14[k]
                  - f_3 * pds1_14[k]
                  + f_4 * pc_x[k] * pdp_42[k];

        t_85[k] = f_0 * ppp_25[k]
                  + f_4 * pc_x[k] * pdp_43[k];

        t_86[k] = f_0 * ppp_26[k]
                  + f_4 * pc_x[k] * pdp_44[k];

        t_87[k] = pb_x[k] * ppd0_51[k]
                  - f_5 * pc_x[k] * ppd1_51[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_z, pb_x, pc_x, pc_y, pc_z, sdd0_18, \
                         sdd1_18, ppd0_53, ppd1_53, pdp_44, pdp_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_4 * pc_y[k] * pdp_44[k];

        t_89[k] = pb_x[k] * ppd0_53[k]
                  - f_5 * pc_x[k] * ppd1_53[k];

        t_90[k] = pa_z[k] * sdd0_18[k]
                  - f_5 * pc_z[k] * sdd1_18[k];

        t_91[k] = f_4 * pc_x[k] * pdp_46[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_z, pc_x, pc_y, pc_z, sdd0_21, sdd0_23, \
                         sdp_11, sdd1_21, sdd1_23, ppp_23, pdp_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_4 * pc_x[k] * pdp_47[k];

        t_93[k] = pa_z[k] * sdd0_21[k]
                  - f_5 * pc_z[k] * sdd1_21[k];

        t_94[k] = f_1 * ppp_23[k]
                  + f_4 * pc_y[k] * pdp_47[k];

        t_95[k] = pa_z[k] * sdd0_23[k]
                  + f_1 * sdp_11[k]
                  - f_5 * pc_z[k] * sdd1_23[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, pb_y, pc_x, pc_y, ppd0_48, ppp_25, \
                         ppp_26, ppd1_48, pds0_16, pds1_16, pdp_49, \
                         pdp_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pb_y[k] * ppd0_48[k]
                  - f_5 * pc_y[k] * ppd1_48[k];

        t_97[k] = f_4 * pc_x[k] * pdp_49[k];

        t_98[k] = f_4 * pc_x[k] * pdp_50[k];

        t_99[k] = f_0 * ppp_25[k]
                  + f_2 * pds0_16[k]
                  - f_3 * pds1_16[k]
                  + f_4 * pc_y[k] * pdp_49[k];

        t_100[k] = f_0 * ppp_26[k]
                   + f_4 * pc_y[k] * pdp_50[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, t_105, t_106, pb_y, pc_x, pc_y, ppd0_53, \
                         ppd1_53, pds0_17, pds1_17, pdp_51, pdp_52, \
                         pdp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = pb_y[k] * ppd0_53[k]
                   - f_5 * pc_y[k] * ppd1_53[k];

        t_102[k] = f_2 * pds0_17[k]
                   - f_3 * pds1_17[k]
                   + f_4 * pc_x[k] * pdp_51[k];

        t_103[k] = f_4 * pc_x[k] * pdp_52[k];

        t_104[k] = f_4 * pc_x[k] * pdp_53[k];

        t_105[k] = f_2 * pds0_17[k]
                   - f_3 * pds1_17[k]
                   + f_4 * pc_y[k] * pdp_52[k];

        t_106[k] = f_4 * pc_y[k] * pdp_53[k];
    }

#pragma omp simd aligned(t_107, pc_z, sdp_17, ppp_26, pds0_17, pds1_17, \
                         pdp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_0 * sdp_17[k]
                   + f_1 * ppp_26[k]
                   + f_2 * pds0_17[k]
                   - f_3 * pds1_17[k]
                   + f_4 * pc_z[k] * pdp_53[k];
    }
}

}  // namespace simdt3ceri
