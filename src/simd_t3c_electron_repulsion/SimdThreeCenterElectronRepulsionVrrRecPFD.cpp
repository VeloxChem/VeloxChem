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


#include "SimdThreeCenterElectronRepulsionVrrRecPFD.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_pfd_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sfd0, const size_t sfp,
                                                          const size_t sfd1, const size_t ppd0,
                                                          const size_t ppd1, const size_t pdd0,
                                                          const size_t pdp, const size_t pdd1,
                                                          const size_t pfs0, const size_t pfs1,
                                                          const size_t pfp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.5 / q;
    const auto f_2 = 0.5 / gamma;
    const auto f_3 = 0.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = gamma / q;
    const auto f_6 = 1.0 / q;
    const auto f_7 = 0.5 / p;
    const auto f_8 = 0.5 * gamma / (p * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfd0_0 = buffer.data(sfd0 + 0);
    const auto *sfd0_3 = buffer.data(sfd0 + 3);
    const auto *sfd0_5 = buffer.data(sfd0 + 5);
    const auto *sfd0_12 = buffer.data(sfd0 + 12);
    const auto *sfd0_17 = buffer.data(sfd0 + 17);
    const auto *sfd0_30 = buffer.data(sfd0 + 30);
    const auto *sfd0_35 = buffer.data(sfd0 + 35);
    const auto *sfd0_36 = buffer.data(sfd0 + 36);
    const auto *sfd0_39 = buffer.data(sfd0 + 39);
    const auto *sfd0_41 = buffer.data(sfd0 + 41);
    const auto *sfd0_45 = buffer.data(sfd0 + 45);
    const auto *sfd0_47 = buffer.data(sfd0 + 47);
    const auto *sfd0_51 = buffer.data(sfd0 + 51);
    const auto *sfd0_53 = buffer.data(sfd0 + 53);
    const auto *sfd0_54 = buffer.data(sfd0 + 54);
    const auto *sfd0_57 = buffer.data(sfd0 + 57);
    const auto *sfd0_59 = buffer.data(sfd0 + 59);

    const auto *sfp_0 = buffer.data(sfp + 0);
    const auto *sfp_1 = buffer.data(sfp + 1);
    const auto *sfp_9 = buffer.data(sfp + 9);
    const auto *sfp_15 = buffer.data(sfp + 15);
    const auto *sfp_18 = buffer.data(sfp + 18);
    const auto *sfp_19 = buffer.data(sfp + 19);
    const auto *sfp_25 = buffer.data(sfp + 25);
    const auto *sfp_27 = buffer.data(sfp + 27);
    const auto *sfp_28 = buffer.data(sfp + 28);

    const auto *sfd1_0 = buffer.data(sfd1 + 0);
    const auto *sfd1_3 = buffer.data(sfd1 + 3);
    const auto *sfd1_5 = buffer.data(sfd1 + 5);
    const auto *sfd1_12 = buffer.data(sfd1 + 12);
    const auto *sfd1_17 = buffer.data(sfd1 + 17);
    const auto *sfd1_30 = buffer.data(sfd1 + 30);
    const auto *sfd1_35 = buffer.data(sfd1 + 35);
    const auto *sfd1_36 = buffer.data(sfd1 + 36);
    const auto *sfd1_39 = buffer.data(sfd1 + 39);
    const auto *sfd1_41 = buffer.data(sfd1 + 41);
    const auto *sfd1_45 = buffer.data(sfd1 + 45);
    const auto *sfd1_47 = buffer.data(sfd1 + 47);
    const auto *sfd1_51 = buffer.data(sfd1 + 51);
    const auto *sfd1_53 = buffer.data(sfd1 + 53);
    const auto *sfd1_54 = buffer.data(sfd1 + 54);
    const auto *sfd1_57 = buffer.data(sfd1 + 57);
    const auto *sfd1_59 = buffer.data(sfd1 + 59);

    const auto *ppd0_27 = buffer.data(ppd0 + 27);

    const auto *ppd1_27 = buffer.data(ppd1 + 27);

    const auto *pdd0_0 = buffer.data(pdd0 + 0);
    const auto *pdd0_3 = buffer.data(pdd0 + 3);
    const auto *pdd0_5 = buffer.data(pdd0 + 5);
    const auto *pdd0_12 = buffer.data(pdd0 + 12);
    const auto *pdd0_18 = buffer.data(pdd0 + 18);
    const auto *pdd0_30 = buffer.data(pdd0 + 30);
    const auto *pdd0_39 = buffer.data(pdd0 + 39);
    const auto *pdd0_45 = buffer.data(pdd0 + 45);
    const auto *pdd0_54 = buffer.data(pdd0 + 54);
    const auto *pdd0_57 = buffer.data(pdd0 + 57);
    const auto *pdd0_59 = buffer.data(pdd0 + 59);
    const auto *pdd0_63 = buffer.data(pdd0 + 63);
    const auto *pdd0_65 = buffer.data(pdd0 + 65);
    const auto *pdd0_69 = buffer.data(pdd0 + 69);

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

    const auto *pdd1_0 = buffer.data(pdd1 + 0);
    const auto *pdd1_3 = buffer.data(pdd1 + 3);
    const auto *pdd1_5 = buffer.data(pdd1 + 5);
    const auto *pdd1_12 = buffer.data(pdd1 + 12);
    const auto *pdd1_18 = buffer.data(pdd1 + 18);
    const auto *pdd1_30 = buffer.data(pdd1 + 30);
    const auto *pdd1_39 = buffer.data(pdd1 + 39);
    const auto *pdd1_45 = buffer.data(pdd1 + 45);
    const auto *pdd1_54 = buffer.data(pdd1 + 54);
    const auto *pdd1_57 = buffer.data(pdd1 + 57);
    const auto *pdd1_59 = buffer.data(pdd1 + 59);
    const auto *pdd1_63 = buffer.data(pdd1 + 63);
    const auto *pdd1_65 = buffer.data(pdd1 + 65);
    const auto *pdd1_69 = buffer.data(pdd1 + 69);

    const auto *pfs0_0 = buffer.data(pfs0 + 0);
    const auto *pfs0_1 = buffer.data(pfs0 + 1);
    const auto *pfs0_2 = buffer.data(pfs0 + 2);
    const auto *pfs0_3 = buffer.data(pfs0 + 3);
    const auto *pfs0_4 = buffer.data(pfs0 + 4);
    const auto *pfs0_5 = buffer.data(pfs0 + 5);
    const auto *pfs0_11 = buffer.data(pfs0 + 11);
    const auto *pfs0_14 = buffer.data(pfs0 + 14);
    const auto *pfs0_16 = buffer.data(pfs0 + 16);
    const auto *pfs0_17 = buffer.data(pfs0 + 17);
    const auto *pfs0_18 = buffer.data(pfs0 + 18);

    const auto *pfs1_0 = buffer.data(pfs1 + 0);
    const auto *pfs1_1 = buffer.data(pfs1 + 1);
    const auto *pfs1_2 = buffer.data(pfs1 + 2);
    const auto *pfs1_3 = buffer.data(pfs1 + 3);
    const auto *pfs1_4 = buffer.data(pfs1 + 4);
    const auto *pfs1_5 = buffer.data(pfs1 + 5);
    const auto *pfs1_11 = buffer.data(pfs1 + 11);
    const auto *pfs1_14 = buffer.data(pfs1 + 14);
    const auto *pfs1_16 = buffer.data(pfs1 + 16);
    const auto *pfs1_17 = buffer.data(pfs1 + 17);
    const auto *pfs1_18 = buffer.data(pfs1 + 18);

    const auto *pfp_0 = buffer.data(pfp + 0);
    const auto *pfp_1 = buffer.data(pfp + 1);
    const auto *pfp_2 = buffer.data(pfp + 2);
    const auto *pfp_3 = buffer.data(pfp + 3);
    const auto *pfp_4 = buffer.data(pfp + 4);
    const auto *pfp_5 = buffer.data(pfp + 5);
    const auto *pfp_6 = buffer.data(pfp + 6);
    const auto *pfp_8 = buffer.data(pfp + 8);
    const auto *pfp_9 = buffer.data(pfp + 9);
    const auto *pfp_10 = buffer.data(pfp + 10);
    const auto *pfp_11 = buffer.data(pfp + 11);
    const auto *pfp_12 = buffer.data(pfp + 12);
    const auto *pfp_13 = buffer.data(pfp + 13);
    const auto *pfp_14 = buffer.data(pfp + 14);
    const auto *pfp_15 = buffer.data(pfp + 15);
    const auto *pfp_16 = buffer.data(pfp + 16);
    const auto *pfp_17 = buffer.data(pfp + 17);
    const auto *pfp_18 = buffer.data(pfp + 18);
    const auto *pfp_20 = buffer.data(pfp + 20);
    const auto *pfp_21 = buffer.data(pfp + 21);
    const auto *pfp_23 = buffer.data(pfp + 23);
    const auto *pfp_24 = buffer.data(pfp + 24);
    const auto *pfp_26 = buffer.data(pfp + 26);
    const auto *pfp_27 = buffer.data(pfp + 27);
    const auto *pfp_29 = buffer.data(pfp + 29);
    const auto *pfp_31 = buffer.data(pfp + 31);
    const auto *pfp_32 = buffer.data(pfp + 32);
    const auto *pfp_33 = buffer.data(pfp + 33);
    const auto *pfp_34 = buffer.data(pfp + 34);
    const auto *pfp_35 = buffer.data(pfp + 35);
    const auto *pfp_37 = buffer.data(pfp + 37);
    const auto *pfp_38 = buffer.data(pfp + 38);
    const auto *pfp_40 = buffer.data(pfp + 40);
    const auto *pfp_41 = buffer.data(pfp + 41);
    const auto *pfp_42 = buffer.data(pfp + 42);
    const auto *pfp_43 = buffer.data(pfp + 43);
    const auto *pfp_44 = buffer.data(pfp + 44);
    const auto *pfp_46 = buffer.data(pfp + 46);
    const auto *pfp_47 = buffer.data(pfp + 47);
    const auto *pfp_48 = buffer.data(pfp + 48);
    const auto *pfp_49 = buffer.data(pfp + 49);
    const auto *pfp_50 = buffer.data(pfp + 50);
    const auto *pfp_52 = buffer.data(pfp + 52);
    const auto *pfp_53 = buffer.data(pfp + 53);
    const auto *pfp_54 = buffer.data(pfp + 54);
    const auto *pfp_55 = buffer.data(pfp + 55);
    const auto *pfp_56 = buffer.data(pfp + 56);
    const auto *pfp_58 = buffer.data(pfp + 58);
    const auto *pfp_59 = buffer.data(pfp + 59);
    const auto *pfp_61 = buffer.data(pfp + 61);
    const auto *pfp_62 = buffer.data(pfp + 62);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, sfp_0, pdp_0, pfs0_0, \
                         pfs1_0, pfp_0, pfp_1, pfp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sfp_0[k]
                 + f_1 * pdp_0[k]
                 + f_2 * pfs0_0[k]
                 - f_3 * pfs1_0[k]
                 + f_4 * pc_x[k] * pfp_0[k];

        t_1[k] = f_4 * pc_y[k] * pfp_0[k];

        t_2[k] = f_4 * pc_z[k] * pfp_0[k];

        t_3[k] = f_2 * pfs0_0[k]
                 - f_3 * pfs1_0[k]
                 + f_4 * pc_y[k] * pfp_1[k];

        t_4[k] = f_4 * pc_y[k] * pfp_2[k];

        t_5[k] = f_2 * pfs0_0[k]
                 - f_3 * pfs1_0[k]
                 + f_4 * pc_z[k] * pfp_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pb_y, pc_y, pc_z, pdd0_0, pdp_0, pdp_1, pdd1_0, \
                         pfs0_1, pfs1_1, pfp_3, pfp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pb_y[k] * pdd0_0[k]
                 - f_5 * pc_y[k] * pdd1_0[k];

        t_7[k] = f_0 * pdp_0[k]
                 + f_4 * pc_y[k] * pfp_3[k];

        t_8[k] = f_4 * pc_z[k] * pfp_3[k];

        t_9[k] = f_0 * pdp_1[k]
                 + f_2 * pfs0_1[k]
                 - f_3 * pfs1_1[k]
                 + f_4 * pc_y[k] * pfp_4[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pb_y, pb_z, pc_y, pc_z, pdd0_0, pdd0_5, \
                         pdp_2, pdd1_0, pdd1_5, pfp_5, pfp_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_0 * pdp_2[k]
                  + f_4 * pc_y[k] * pfp_5[k];

        t_11[k] = pb_y[k] * pdd0_5[k]
                  - f_5 * pc_y[k] * pdd1_5[k];

        t_12[k] = pb_z[k] * pdd0_0[k]
                  - f_5 * pc_z[k] * pdd1_0[k];

        t_13[k] = f_4 * pc_y[k] * pfp_6[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pb_z, pc_y, pc_z, pdd0_3, pdp_0, pdp_2, \
                         pdd1_3, pfs0_2, pfs1_2, pfp_6, pfp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_0 * pdp_0[k]
                  + f_4 * pc_z[k] * pfp_6[k];

        t_15[k] = pb_z[k] * pdd0_3[k]
                  - f_5 * pc_z[k] * pdd1_3[k];

        t_16[k] = f_4 * pc_y[k] * pfp_8[k];

        t_17[k] = f_0 * pdp_2[k]
                  + f_2 * pfs0_2[k]
                  - f_3 * pfs1_2[k]
                  + f_4 * pc_z[k] * pfp_8[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pc_x, pc_y, pc_z, sfp_9, pdp_3, pdp_4, pdp_9, \
                         pfs0_3, pfs1_3, pfp_9, pfp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_0 * sfp_9[k]
                  + f_0 * pdp_9[k]
                  + f_2 * pfs0_3[k]
                  - f_3 * pfs1_3[k]
                  + f_4 * pc_x[k] * pfp_9[k];

        t_19[k] = f_6 * pdp_3[k]
                  + f_4 * pc_y[k] * pfp_9[k];

        t_20[k] = f_4 * pc_z[k] * pfp_9[k];

        t_21[k] = f_6 * pdp_4[k]
                  + f_2 * pfs0_3[k]
                  - f_3 * pfs1_3[k]
                  + f_4 * pc_y[k] * pfp_10[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pb_y, pc_y, pc_z, pdd0_12, pdp_5, pdp_6, \
                         pdd1_12, pfs0_3, pfs1_3, pfp_11, pfp_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_6 * pdp_5[k]
                  + f_4 * pc_y[k] * pfp_11[k];

        t_23[k] = f_2 * pfs0_3[k]
                  - f_3 * pfs1_3[k]
                  + f_4 * pc_z[k] * pfp_11[k];

        t_24[k] = pb_y[k] * pdd0_12[k]
                  - f_5 * pc_y[k] * pdd1_12[k];

        t_25[k] = f_0 * pdp_6[k]
                  + f_4 * pc_y[k] * pfp_12[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pc_y, pc_z, pdp_3, pdp_5, pdp_7, pdp_8, \
                         pfs0_4, pfs1_4, pfp_12, pfp_13, pfp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_0 * pdp_3[k]
                  + f_4 * pc_z[k] * pfp_12[k];

        t_27[k] = f_0 * pdp_7[k]
                  + f_2 * pfs0_4[k]
                  - f_3 * pfs1_4[k]
                  + f_4 * pc_y[k] * pfp_13[k];

        t_28[k] = f_0 * pdp_8[k]
                  + f_4 * pc_y[k] * pfp_14[k];

        t_29[k] = f_0 * pdp_5[k]
                  + f_2 * pfs0_4[k]
                  - f_3 * pfs1_4[k]
                  + f_4 * pc_z[k] * pfp_14[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pc_x, pc_y, pc_z, sfp_15, pdp_6, \
                         pdp_15, pfs0_5, pfs1_5, pfp_15, pfp_16, \
                         pfp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_0 * sfp_15[k]
                  + f_0 * pdp_15[k]
                  + f_2 * pfs0_5[k]
                  - f_3 * pfs1_5[k]
                  + f_4 * pc_x[k] * pfp_15[k];

        t_31[k] = f_4 * pc_y[k] * pfp_15[k];

        t_32[k] = f_6 * pdp_6[k]
                  + f_4 * pc_z[k] * pfp_15[k];

        t_33[k] = f_2 * pfs0_5[k]
                  - f_3 * pfs1_5[k]
                  + f_4 * pc_y[k] * pfp_16[k];

        t_34[k] = f_4 * pc_y[k] * pfp_17[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pa_x, pc_x, pc_y, pc_z, sfd0_36, sfp_18, sfd1_36, \
                         pdp_8, pdp_9, pfs0_5, pfs1_5, pfp_17, pfp_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_6 * pdp_8[k]
                  + f_2 * pfs0_5[k]
                  - f_3 * pfs1_5[k]
                  + f_4 * pc_z[k] * pfp_17[k];

        t_36[k] = pa_x[k] * sfd0_36[k]
                  + f_6 * sfp_18[k]
                  - f_5 * pc_x[k] * sfd1_36[k];

        t_37[k] = f_1 * pdp_9[k]
                  + f_4 * pc_y[k] * pfp_18[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pa_x, pc_x, pc_y, pc_z, sfd0_39, sfd0_41, \
                         sfd1_39, sfd1_41, pdp_11, pfp_18, pfp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_4 * pc_z[k] * pfp_18[k];

        t_39[k] = pa_x[k] * sfd0_39[k]
                  - f_5 * pc_x[k] * sfd1_39[k];

        t_40[k] = f_1 * pdp_11[k]
                  + f_4 * pc_y[k] * pfp_20[k];

        t_41[k] = pa_x[k] * sfd0_41[k]
                  - f_5 * pc_x[k] * sfd1_41[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_x, pb_z, pc_x, pc_y, pc_z, sfd0_45, \
                         sfd1_45, pdd0_18, pdp_9, pdp_12, pdd1_18, \
                         pfp_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_z[k] * pdd0_18[k]
                  - f_5 * pc_z[k] * pdd1_18[k];

        t_43[k] = f_6 * pdp_12[k]
                  + f_4 * pc_y[k] * pfp_21[k];

        t_44[k] = f_0 * pdp_9[k]
                  + f_4 * pc_z[k] * pfp_21[k];

        t_45[k] = pa_x[k] * sfd0_45[k]
                  - f_5 * pc_x[k] * sfd1_45[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_x, pb_y, pc_x, pc_y, sfd0_47, sfd1_47, \
                         pdd0_30, pdp_14, pdp_15, pdd1_30, pfp_23, \
                         pfp_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_6 * pdp_14[k]
                  + f_4 * pc_y[k] * pfp_23[k];

        t_47[k] = pa_x[k] * sfd0_47[k]
                  - f_5 * pc_x[k] * sfd1_47[k];

        t_48[k] = pb_y[k] * pdd0_30[k]
                  - f_5 * pc_y[k] * pdd1_30[k];

        t_49[k] = f_0 * pdp_15[k]
                  + f_4 * pc_y[k] * pfp_24[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_x, pc_x, pc_y, pc_z, sfd0_51, sfd0_53, \
                         sfd1_51, sfd1_53, pdp_12, pdp_17, pfp_24, \
                         pfp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_6 * pdp_12[k]
                  + f_4 * pc_z[k] * pfp_24[k];

        t_51[k] = pa_x[k] * sfd0_51[k]
                  - f_5 * pc_x[k] * sfd1_51[k];

        t_52[k] = f_0 * pdp_17[k]
                  + f_4 * pc_y[k] * pfp_26[k];

        t_53[k] = pa_x[k] * sfd0_53[k]
                  - f_5 * pc_x[k] * sfd1_53[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_x, pc_x, pc_y, pc_z, sfd0_54, sfd0_57, \
                         sfp_27, sfd1_54, sfd1_57, pdp_15, pfp_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = pa_x[k] * sfd0_54[k]
                  + f_6 * sfp_27[k]
                  - f_5 * pc_x[k] * sfd1_54[k];

        t_55[k] = f_4 * pc_y[k] * pfp_27[k];

        t_56[k] = f_1 * pdp_15[k]
                  + f_4 * pc_z[k] * pfp_27[k];

        t_57[k] = pa_x[k] * sfd0_57[k]
                  - f_5 * pc_x[k] * sfd1_57[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pa_x, pa_y, pc_x, pc_y, sfd0_0, sfd0_59, \
                         sfd1_0, sfd1_59, pdp_19, pfp_29, pfp_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_4 * pc_y[k] * pfp_29[k];

        t_59[k] = pa_x[k] * sfd0_59[k]
                  - f_5 * pc_x[k] * sfd1_59[k];

        t_60[k] = pa_y[k] * sfd0_0[k]
                  - f_5 * pc_y[k] * sfd1_0[k];

        t_61[k] = f_1 * pdp_19[k]
                  + f_4 * pc_x[k] * pfp_31[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pa_y, pc_x, pc_y, pc_z, sfd0_3, sfd0_5, \
                         sfp_1, sfd1_3, sfd1_5, pdp_20, pfp_31, \
                         pfp_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_1 * pdp_20[k]
                  + f_4 * pc_x[k] * pfp_32[k];

        t_63[k] = pa_y[k] * sfd0_3[k]
                  + f_6 * sfp_1[k]
                  - f_5 * pc_y[k] * sfd1_3[k];

        t_64[k] = f_4 * pc_z[k] * pfp_31[k];

        t_65[k] = pa_y[k] * sfd0_5[k]
                  - f_5 * pc_y[k] * sfd1_5[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pc_x, pdp_21, pdp_22, pdp_23, pfs0_11, pfs1_11, \
                         pfp_33, pfp_34, pfp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_6 * pdp_21[k]
                  + f_2 * pfs0_11[k]
                  - f_3 * pfs1_11[k]
                  + f_4 * pc_x[k] * pfp_33[k];

        t_67[k] = f_6 * pdp_22[k]
                  + f_4 * pc_x[k] * pfp_34[k];

        t_68[k] = f_6 * pdp_23[k]
                  + f_4 * pc_x[k] * pfp_35[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pb_x, pc_x, pc_z, ppd0_27, ppd1_27, pdd0_45, \
                         pdd1_45, pfs0_11, pfs1_11, pfp_34, pfp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_7 * ppd0_27[k]
                  - f_8 * ppd1_27[k]
                  + pb_x[k] * pdd0_45[k]
                  - f_5 * pc_x[k] * pdd1_45[k];

        t_70[k] = f_4 * pc_z[k] * pfp_34[k];

        t_71[k] = f_2 * pfs0_11[k]
                  - f_3 * pfs1_11[k]
                  + f_4 * pc_z[k] * pfp_35[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pa_y, pc_x, pc_y, sfd0_12, sfd1_12, pdp_25, pdp_26, \
                         pfp_37, pfp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = pa_y[k] * sfd0_12[k]
                  - f_5 * pc_y[k] * sfd1_12[k];

        t_73[k] = f_6 * pdp_25[k]
                  + f_4 * pc_x[k] * pfp_37[k];

        t_74[k] = f_6 * pdp_26[k]
                  + f_4 * pc_x[k] * pfp_38[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, pa_y, pb_z, pc_y, pc_z, sfd0_17, sfd1_17, pdd0_39, \
                         pdp_19, pdd1_39, pfp_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = pb_z[k] * pdd0_39[k]
                  - f_5 * pc_z[k] * pdd1_39[k];

        t_76[k] = f_0 * pdp_19[k]
                  + f_4 * pc_z[k] * pfp_37[k];

        t_77[k] = pa_y[k] * sfd0_17[k]
                  - f_5 * pc_y[k] * sfd1_17[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pb_x, pc_x, pdd0_54, pdd0_57, pdp_27, pdp_28, \
                         pdp_29, pdd1_54, pdd1_57, pfp_40, pfp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = pb_x[k] * pdd0_54[k]
                  + f_6 * pdp_27[k]
                  - f_5 * pc_x[k] * pdd1_54[k];

        t_79[k] = f_0 * pdp_28[k]
                  + f_4 * pc_x[k] * pfp_40[k];

        t_80[k] = f_0 * pdp_29[k]
                  + f_4 * pc_x[k] * pfp_41[k];

        t_81[k] = pb_x[k] * pdd0_57[k]
                  - f_5 * pc_x[k] * pdd1_57[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pb_x, pc_x, pc_z, pdd0_59, pdp_30, pdp_31, \
                         pdd1_59, pfs0_14, pfs1_14, pfp_40, pfp_42, \
                         pfp_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_4 * pc_z[k] * pfp_40[k];

        t_83[k] = pb_x[k] * pdd0_59[k]
                  - f_5 * pc_x[k] * pdd1_59[k];

        t_84[k] = f_0 * pdp_30[k]
                  + f_2 * pfs0_14[k]
                  - f_3 * pfs1_14[k]
                  + f_4 * pc_x[k] * pfp_42[k];

        t_85[k] = f_0 * pdp_31[k]
                  + f_4 * pc_x[k] * pfp_43[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pb_x, pc_x, pc_z, pdd0_63, pdd0_65, pdp_22, \
                         pdp_32, pdd1_63, pdd1_65, pfp_43, pfp_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_0 * pdp_32[k]
                  + f_4 * pc_x[k] * pfp_44[k];

        t_87[k] = pb_x[k] * pdd0_63[k]
                  - f_5 * pc_x[k] * pdd1_63[k];

        t_88[k] = f_0 * pdp_22[k]
                  + f_4 * pc_z[k] * pfp_43[k];

        t_89[k] = pb_x[k] * pdd0_65[k]
                  - f_5 * pc_x[k] * pdd1_65[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pa_y, pb_x, pc_x, pc_y, sfd0_30, sfd1_30, \
                         pdd0_69, pdp_34, pdp_35, pdd1_69, pfp_46, \
                         pfp_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = pa_y[k] * sfd0_30[k]
                  - f_5 * pc_y[k] * sfd1_30[k];

        t_91[k] = f_0 * pdp_34[k]
                  + f_4 * pc_x[k] * pfp_46[k];

        t_92[k] = f_0 * pdp_35[k]
                  + f_4 * pc_x[k] * pfp_47[k];

        t_93[k] = pb_x[k] * pdd0_69[k]
                  - f_5 * pc_x[k] * pdd1_69[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, pa_y, pc_x, pc_y, pc_z, sfd0_35, sfd1_35, \
                         pdp_25, pfs0_16, pfs1_16, pfp_46, pfp_48, \
                         pfp_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_6 * pdp_25[k]
                  + f_4 * pc_z[k] * pfp_46[k];

        t_95[k] = pa_y[k] * sfd0_35[k]
                  - f_5 * pc_y[k] * sfd1_35[k];

        t_96[k] = f_2 * pfs0_16[k]
                  - f_3 * pfs1_16[k]
                  + f_4 * pc_x[k] * pfp_48[k];

        t_97[k] = f_4 * pc_x[k] * pfp_49[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, t_101, pc_x, pc_y, pc_z, sfp_19, pdp_28, pfs0_16, \
                         pfs1_16, pfp_49, pfp_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_4 * pc_x[k] * pfp_50[k];

        t_99[k] = f_0 * sfp_19[k]
                  + f_1 * pdp_28[k]
                  + f_2 * pfs0_16[k]
                  - f_3 * pfs1_16[k]
                  + f_4 * pc_y[k] * pfp_49[k];

        t_100[k] = f_4 * pc_z[k] * pfp_49[k];

        t_101[k] = f_2 * pfs0_16[k]
                   - f_3 * pfs1_16[k]
                   + f_4 * pc_z[k] * pfp_50[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, pb_z, pc_x, pc_z, pdd0_54, \
                         pdd0_57, pdp_28, pdd1_54, pdd1_57, pfp_52, \
                         pfp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = pb_z[k] * pdd0_54[k]
                   - f_5 * pc_z[k] * pdd1_54[k];

        t_103[k] = f_4 * pc_x[k] * pfp_52[k];

        t_104[k] = f_4 * pc_x[k] * pfp_53[k];

        t_105[k] = pb_z[k] * pdd0_57[k]
                   - f_5 * pc_z[k] * pdd1_57[k];

        t_106[k] = f_0 * pdp_28[k]
                   + f_4 * pc_z[k] * pfp_52[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pc_x, pc_z, pdp_29, pfs0_17, pfs0_18, \
                         pfs1_17, pfs1_18, pfp_53, pfp_54, pfp_55, \
                         pfp_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_0 * pdp_29[k]
                   + f_2 * pfs0_17[k]
                   - f_3 * pfs1_17[k]
                   + f_4 * pc_z[k] * pfp_53[k];

        t_108[k] = f_2 * pfs0_18[k]
                   - f_3 * pfs1_18[k]
                   + f_4 * pc_x[k] * pfp_54[k];

        t_109[k] = f_4 * pc_x[k] * pfp_55[k];

        t_110[k] = f_4 * pc_x[k] * pfp_56[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pc_y, pc_z, sfp_25, pdp_31, pdp_32, pdp_34, \
                         pfs0_18, pfs1_18, pfp_55, pfp_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_0 * sfp_25[k]
                   + f_0 * pdp_34[k]
                   + f_2 * pfs0_18[k]
                   - f_3 * pfs1_18[k]
                   + f_4 * pc_y[k] * pfp_55[k];

        t_112[k] = f_6 * pdp_31[k]
                   + f_4 * pc_z[k] * pfp_55[k];

        t_113[k] = f_6 * pdp_32[k]
                   + f_2 * pfs0_18[k]
                   - f_3 * pfs1_18[k]
                   + f_4 * pc_z[k] * pfp_56[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pa_y, pc_x, pc_y, sfd0_54, sfd0_57, \
                         sfp_28, sfd1_54, sfd1_57, pfp_58, pfp_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = pa_y[k] * sfd0_54[k]
                   - f_5 * pc_y[k] * sfd1_54[k];

        t_115[k] = f_4 * pc_x[k] * pfp_58[k];

        t_116[k] = f_4 * pc_x[k] * pfp_59[k];

        t_117[k] = pa_y[k] * sfd0_57[k]
                   + f_6 * sfp_28[k]
                   - f_5 * pc_y[k] * sfd1_57[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pa_y, pa_z, pc_y, pc_z, sfd0_0, sfd0_59, sfd1_0, \
                         sfd1_59, pdp_34, pfp_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_1 * pdp_34[k]
                   + f_4 * pc_z[k] * pfp_58[k];

        t_119[k] = pa_y[k] * sfd0_59[k]
                   - f_5 * pc_y[k] * sfd1_59[k];

        t_120[k] = pa_z[k] * sfd0_0[k]
                   - f_5 * pc_z[k] * sfd1_0[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pa_z, pc_x, pc_y, pc_z, sfd0_3, sfd1_3, \
                         pdp_37, pdp_38, pfp_61, pfp_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_1 * pdp_37[k]
                   + f_4 * pc_x[k] * pfp_61[k];

        t_122[k] = f_1 * pdp_38[k]
                   + f_4 * pc_x[k] * pfp_62[k];

        t_123[k] = pa_z[k] * sfd0_3[k]
                   - f_5 * pc_z[k] * sfd1_3[k];

        t_124[k] = f_4 * pc_y[k] * pfp_62[k];
    }
}

static auto
compute_prim_pfd_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sfd0, const size_t sfp,
                                                          const size_t sfd1, const size_t ppd0,
                                                          const size_t ppd1, const size_t pdd0,
                                                          const size_t pdp, const size_t pdd1,
                                                          const size_t pfs0, const size_t pfs1,
                                                          const size_t pfp, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.5 / q;
    const auto f_2 = 0.5 / gamma;
    const auto f_3 = 0.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = gamma / q;
    const auto f_6 = 1.0 / q;
    const auto f_7 = 0.5 / p;
    const auto f_8 = 0.5 * gamma / (p * q);

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

    const auto *sfd0_5 = buffer.data(sfd0 + 5);
    const auto *sfd0_6 = buffer.data(sfd0 + 6);
    const auto *sfd0_9 = buffer.data(sfd0 + 9);
    const auto *sfd0_18 = buffer.data(sfd0 + 18);
    const auto *sfd0_21 = buffer.data(sfd0 + 21);
    const auto *sfd0_36 = buffer.data(sfd0 + 36);
    const auto *sfd0_39 = buffer.data(sfd0 + 39);
    const auto *sfd0_41 = buffer.data(sfd0 + 41);

    const auto *sfp_2 = buffer.data(sfp + 2);
    const auto *sfp_20 = buffer.data(sfp + 20);
    const auto *sfp_29 = buffer.data(sfp + 29);

    const auto *sfd1_5 = buffer.data(sfd1 + 5);
    const auto *sfd1_6 = buffer.data(sfd1 + 6);
    const auto *sfd1_9 = buffer.data(sfd1 + 9);
    const auto *sfd1_18 = buffer.data(sfd1 + 18);
    const auto *sfd1_21 = buffer.data(sfd1 + 21);
    const auto *sfd1_36 = buffer.data(sfd1 + 36);
    const auto *sfd1_39 = buffer.data(sfd1 + 39);
    const auto *sfd1_41 = buffer.data(sfd1 + 41);

    const auto *ppd0_53 = buffer.data(ppd0 + 53);

    const auto *ppd1_53 = buffer.data(ppd1 + 53);

    const auto *pdd0_77 = buffer.data(pdd0 + 77);
    const auto *pdd0_84 = buffer.data(pdd0 + 84);
    const auto *pdd0_89 = buffer.data(pdd0 + 89);
    const auto *pdd0_95 = buffer.data(pdd0 + 95);
    const auto *pdd0_99 = buffer.data(pdd0 + 99);
    const auto *pdd0_101 = buffer.data(pdd0 + 101);
    const auto *pdd0_102 = buffer.data(pdd0 + 102);
    const auto *pdd0_105 = buffer.data(pdd0 + 105);
    const auto *pdd0_107 = buffer.data(pdd0 + 107);

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

    const auto *pdd1_77 = buffer.data(pdd1 + 77);
    const auto *pdd1_84 = buffer.data(pdd1 + 84);
    const auto *pdd1_89 = buffer.data(pdd1 + 89);
    const auto *pdd1_95 = buffer.data(pdd1 + 95);
    const auto *pdd1_99 = buffer.data(pdd1 + 99);
    const auto *pdd1_101 = buffer.data(pdd1 + 101);
    const auto *pdd1_102 = buffer.data(pdd1 + 102);
    const auto *pdd1_105 = buffer.data(pdd1 + 105);
    const auto *pdd1_107 = buffer.data(pdd1 + 107);

    const auto *pfs0_22 = buffer.data(pfs0 + 22);
    const auto *pfs0_27 = buffer.data(pfs0 + 27);
    const auto *pfs0_29 = buffer.data(pfs0 + 29);

    const auto *pfs1_22 = buffer.data(pfs1 + 22);
    const auto *pfs1_27 = buffer.data(pfs1 + 27);
    const auto *pfs1_29 = buffer.data(pfs1 + 29);

    const auto *pfp_64 = buffer.data(pfp + 64);
    const auto *pfp_65 = buffer.data(pfp + 65);
    const auto *pfp_66 = buffer.data(pfp + 66);
    const auto *pfp_67 = buffer.data(pfp + 67);
    const auto *pfp_68 = buffer.data(pfp + 68);
    const auto *pfp_70 = buffer.data(pfp + 70);
    const auto *pfp_71 = buffer.data(pfp + 71);
    const auto *pfp_73 = buffer.data(pfp + 73);
    const auto *pfp_74 = buffer.data(pfp + 74);
    const auto *pfp_76 = buffer.data(pfp + 76);
    const auto *pfp_77 = buffer.data(pfp + 77);
    const auto *pfp_79 = buffer.data(pfp + 79);
    const auto *pfp_80 = buffer.data(pfp + 80);
    const auto *pfp_81 = buffer.data(pfp + 81);
    const auto *pfp_82 = buffer.data(pfp + 82);
    const auto *pfp_83 = buffer.data(pfp + 83);
    const auto *pfp_85 = buffer.data(pfp + 85);
    const auto *pfp_86 = buffer.data(pfp + 86);
    const auto *pfp_87 = buffer.data(pfp + 87);
    const auto *pfp_88 = buffer.data(pfp + 88);
    const auto *pfp_89 = buffer.data(pfp + 89);

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pa_z, pc_x, pc_z, sfd0_5, sfd0_6, sfp_2, \
                         sfd1_5, sfd1_6, pdp_40, pdp_41, pfp_64, \
                         pfp_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = pa_z[k] * sfd0_5[k]
                   + f_6 * sfp_2[k]
                   - f_5 * pc_z[k] * sfd1_5[k];

        t_126[k] = pa_z[k] * sfd0_6[k]
                   - f_5 * pc_z[k] * sfd1_6[k];

        t_127[k] = f_6 * pdp_40[k]
                   + f_4 * pc_x[k] * pfp_64[k];

        t_128[k] = f_6 * pdp_41[k]
                   + f_4 * pc_x[k] * pfp_65[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, pa_z, pb_y, pc_y, pc_z, sfd0_9, sfd1_9, pdd0_77, \
                         pdp_38, pdd1_77, pfp_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = pa_z[k] * sfd0_9[k]
                   - f_5 * pc_z[k] * sfd1_9[k];

        t_130[k] = f_0 * pdp_38[k]
                   + f_4 * pc_y[k] * pfp_65[k];

        t_131[k] = pb_y[k] * pdd0_77[k]
                   - f_5 * pc_y[k] * pdd1_77[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, pc_x, pc_y, pdp_42, pdp_43, \
                         pdp_44, pfs0_22, pfs1_22, pfp_66, pfp_67, \
                         pfp_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_6 * pdp_42[k]
                   + f_2 * pfs0_22[k]
                   - f_3 * pfs1_22[k]
                   + f_4 * pc_x[k] * pfp_66[k];

        t_133[k] = f_6 * pdp_43[k]
                   + f_4 * pc_x[k] * pfp_67[k];

        t_134[k] = f_6 * pdp_44[k]
                   + f_4 * pc_x[k] * pfp_68[k];

        t_135[k] = f_2 * pfs0_22[k]
                   - f_3 * pfs1_22[k]
                   + f_4 * pc_y[k] * pfp_67[k];

        t_136[k] = f_4 * pc_y[k] * pfp_68[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pa_z, pb_x, pc_x, pc_z, sfd0_18, sfd1_18, \
                         ppd0_53, ppd1_53, pdd0_89, pdp_46, pdd1_89, \
                         pfp_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_7 * ppd0_53[k]
                   - f_8 * ppd1_53[k]
                   + pb_x[k] * pdd0_89[k]
                   - f_5 * pc_x[k] * pdd1_89[k];

        t_138[k] = pa_z[k] * sfd0_18[k]
                   - f_5 * pc_z[k] * sfd1_18[k];

        t_139[k] = f_0 * pdp_46[k]
                   + f_4 * pc_x[k] * pfp_70[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pa_z, pb_x, pc_x, pc_y, pc_z, sfd0_21, \
                         sfd1_21, pdd0_95, pdp_41, pdp_47, pdd1_95, \
                         pfp_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_0 * pdp_47[k]
                   + f_4 * pc_x[k] * pfp_71[k];

        t_141[k] = pa_z[k] * sfd0_21[k]
                   - f_5 * pc_z[k] * sfd1_21[k];

        t_142[k] = f_6 * pdp_41[k]
                   + f_4 * pc_y[k] * pfp_71[k];

        t_143[k] = pb_x[k] * pdd0_95[k]
                   - f_5 * pc_x[k] * pdd1_95[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pb_x, pb_y, pc_x, pc_y, pdd0_84, pdd0_99, \
                         pdp_49, pdp_50, pdd1_84, pdd1_99, pfp_73, \
                         pfp_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = pb_y[k] * pdd0_84[k]
                   - f_5 * pc_y[k] * pdd1_84[k];

        t_145[k] = f_0 * pdp_49[k]
                   + f_4 * pc_x[k] * pfp_73[k];

        t_146[k] = f_0 * pdp_50[k]
                   + f_4 * pc_x[k] * pfp_74[k];

        t_147[k] = pb_x[k] * pdd0_99[k]
                   - f_5 * pc_x[k] * pdd1_99[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, pb_x, pc_x, pc_y, pdd0_101, pdd0_102, \
                         pdp_44, pdp_51, pdp_52, pdd1_101, pdd1_102, pfp_74, \
                         pfp_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_0 * pdp_44[k]
                   + f_4 * pc_y[k] * pfp_74[k];

        t_149[k] = pb_x[k] * pdd0_101[k]
                   - f_5 * pc_x[k] * pdd1_101[k];

        t_150[k] = pb_x[k] * pdd0_102[k]
                   + f_6 * pdp_51[k]
                   - f_5 * pc_x[k] * pdd1_102[k];

        t_151[k] = f_0 * pdp_52[k]
                   + f_4 * pc_x[k] * pfp_76[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pb_x, pc_x, pc_y, pdd0_105, pdd0_107, \
                         pdp_53, pdd1_105, pdd1_107, pfp_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_0 * pdp_53[k]
                   + f_4 * pc_x[k] * pfp_77[k];

        t_153[k] = pb_x[k] * pdd0_105[k]
                   - f_5 * pc_x[k] * pdd1_105[k];

        t_154[k] = f_4 * pc_y[k] * pfp_77[k];

        t_155[k] = pb_x[k] * pdd0_107[k]
                   - f_5 * pc_x[k] * pdd1_107[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, t_160, pa_z, pc_x, pc_y, pc_z, sfd0_36, \
                         sfd0_39, sfd1_36, sfd1_39, pdp_47, pfp_79, \
                         pfp_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = pa_z[k] * sfd0_36[k]
                   - f_5 * pc_z[k] * sfd1_36[k];

        t_157[k] = f_4 * pc_x[k] * pfp_79[k];

        t_158[k] = f_4 * pc_x[k] * pfp_80[k];

        t_159[k] = pa_z[k] * sfd0_39[k]
                   - f_5 * pc_z[k] * sfd1_39[k];

        t_160[k] = f_1 * pdp_47[k]
                   + f_4 * pc_y[k] * pfp_80[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, pa_z, pc_x, pc_z, sfd0_41, sfp_20, \
                         sfd1_41, pfs0_27, pfs1_27, pfp_81, pfp_82, \
                         pfp_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = pa_z[k] * sfd0_41[k]
                   + f_6 * sfp_20[k]
                   - f_5 * pc_z[k] * sfd1_41[k];

        t_162[k] = f_2 * pfs0_27[k]
                   - f_3 * pfs1_27[k]
                   + f_4 * pc_x[k] * pfp_81[k];

        t_163[k] = f_4 * pc_x[k] * pfp_82[k];

        t_164[k] = f_4 * pc_x[k] * pfp_83[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, pb_y, pc_y, ppd0_53, ppd1_53, pdd0_101, pdp_49, \
                         pdp_50, pdd1_101, pfs0_27, pfs1_27, pfp_82, \
                         pfp_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_6 * pdp_49[k]
                   + f_2 * pfs0_27[k]
                   - f_3 * pfs1_27[k]
                   + f_4 * pc_y[k] * pfp_82[k];

        t_166[k] = f_6 * pdp_50[k]
                   + f_4 * pc_y[k] * pfp_83[k];

        t_167[k] = f_7 * ppd0_53[k]
                   - f_8 * ppd1_53[k]
                   + pb_y[k] * pdd0_101[k]
                   - f_5 * pc_y[k] * pdd1_101[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, t_172, pb_y, pc_x, pc_y, pdd0_102, \
                         pdd0_105, pdp_52, pdp_53, pdd1_102, pdd1_105, pfp_85, \
                         pfp_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = pb_y[k] * pdd0_102[k]
                   - f_5 * pc_y[k] * pdd1_102[k];

        t_169[k] = f_4 * pc_x[k] * pfp_85[k];

        t_170[k] = f_4 * pc_x[k] * pfp_86[k];

        t_171[k] = pb_y[k] * pdd0_105[k]
                   + f_6 * pdp_52[k]
                   - f_5 * pc_y[k] * pdd1_105[k];

        t_172[k] = f_0 * pdp_53[k]
                   + f_4 * pc_y[k] * pfp_86[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, t_177, t_178, pb_y, pc_x, pc_y, pdd0_107, \
                         pdd1_107, pfs0_29, pfs1_29, pfp_87, pfp_88, \
                         pfp_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = pb_y[k] * pdd0_107[k]
                   - f_5 * pc_y[k] * pdd1_107[k];

        t_174[k] = f_2 * pfs0_29[k]
                   - f_3 * pfs1_29[k]
                   + f_4 * pc_x[k] * pfp_87[k];

        t_175[k] = f_4 * pc_x[k] * pfp_88[k];

        t_176[k] = f_4 * pc_x[k] * pfp_89[k];

        t_177[k] = f_2 * pfs0_29[k]
                   - f_3 * pfs1_29[k]
                   + f_4 * pc_y[k] * pfp_88[k];

        t_178[k] = f_4 * pc_y[k] * pfp_89[k];
    }

#pragma omp simd aligned(t_179, pc_z, sfp_29, pdp_53, pfs0_29, pfs1_29, \
                         pfp_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_0 * sfp_29[k]
                   + f_1 * pdp_53[k]
                   + f_2 * pfs0_29[k]
                   - f_3 * pfs1_29[k]
                   + f_4 * pc_z[k] * pfp_89[k];
    }
}

auto
compute_prim_pfd_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t sfd0,
                                                   const size_t sfp, const size_t sfd1,
                                                   const size_t ppd0, const size_t ppd1,
                                                   const size_t pdd0, const size_t pdp,
                                                   const size_t pdd1, const size_t pfs0,
                                                   const size_t pfs1, const size_t pfp,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_pfd_three_center_electron_repulsion_0_piece0(buffer, target, pa, pb, pc, sfd0,
                                                              sfp, sfd1, ppd0, ppd1, pdd0, pdp,
                                                              pdd1, pfs0, pfs1, pfp, ncols,
                                                              gamma, p, q);

    compute_prim_pfd_three_center_electron_repulsion_0_piece1(buffer, target, pa, pb, pc, sfd0,
                                                              sfp, sfd1, ppd0, ppd1, pdd0, pdp,
                                                              pdd1, pfs0, pfs1, pfp, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
