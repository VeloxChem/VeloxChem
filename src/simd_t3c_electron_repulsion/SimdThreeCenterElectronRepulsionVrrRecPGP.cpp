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


#include "SimdThreeCenterElectronRepulsionVrrRecPGP.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_pgp_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sgs, const size_t pfp0,
                                                   const size_t pfs, const size_t pfp1,
                                                   const size_t pgs, const size_t ncols,
                                                   const double gamma, const double p,
                                                   const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 2.0 / q;
    const auto f_2 = p / q;
    const auto f_3 = gamma / q;
    const auto f_4 = 1.0 / q;
    const auto f_5 = 1.5 / q;

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
    auto *t_133 = buffer.data(target + 133);
    auto *t_134 = buffer.data(target + 134);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sgs_0 = buffer.data(sgs + 0);
    const auto *sgs_1 = buffer.data(sgs + 1);
    const auto *sgs_2 = buffer.data(sgs + 2);
    const auto *sgs_3 = buffer.data(sgs + 3);
    const auto *sgs_5 = buffer.data(sgs + 5);
    const auto *sgs_6 = buffer.data(sgs + 6);
    const auto *sgs_9 = buffer.data(sgs + 9);
    const auto *sgs_10 = buffer.data(sgs + 10);
    const auto *sgs_11 = buffer.data(sgs + 11);
    const auto *sgs_12 = buffer.data(sgs + 12);
    const auto *sgs_13 = buffer.data(sgs + 13);
    const auto *sgs_14 = buffer.data(sgs + 14);

    const auto *pfp0_0 = buffer.data(pfp0 + 0);
    const auto *pfp0_6 = buffer.data(pfp0 + 6);
    const auto *pfp0_9 = buffer.data(pfp0 + 9);
    const auto *pfp0_15 = buffer.data(pfp0 + 15);
    const auto *pfp0_34 = buffer.data(pfp0 + 34);
    const auto *pfp0_49 = buffer.data(pfp0 + 49);
    const auto *pfp0_52 = buffer.data(pfp0 + 52);
    const auto *pfp0_55 = buffer.data(pfp0 + 55);
    const auto *pfp0_68 = buffer.data(pfp0 + 68);
    const auto *pfp0_83 = buffer.data(pfp0 + 83);
    const auto *pfp0_86 = buffer.data(pfp0 + 86);
    const auto *pfp0_89 = buffer.data(pfp0 + 89);

    const auto *pfs_0 = buffer.data(pfs + 0);
    const auto *pfs_1 = buffer.data(pfs + 1);
    const auto *pfs_2 = buffer.data(pfs + 2);
    const auto *pfs_3 = buffer.data(pfs + 3);
    const auto *pfs_4 = buffer.data(pfs + 4);
    const auto *pfs_5 = buffer.data(pfs + 5);
    const auto *pfs_6 = buffer.data(pfs + 6);
    const auto *pfs_7 = buffer.data(pfs + 7);
    const auto *pfs_8 = buffer.data(pfs + 8);
    const auto *pfs_9 = buffer.data(pfs + 9);
    const auto *pfs_10 = buffer.data(pfs + 10);
    const auto *pfs_11 = buffer.data(pfs + 11);
    const auto *pfs_12 = buffer.data(pfs + 12);
    const auto *pfs_13 = buffer.data(pfs + 13);
    const auto *pfs_14 = buffer.data(pfs + 14);
    const auto *pfs_15 = buffer.data(pfs + 15);
    const auto *pfs_16 = buffer.data(pfs + 16);
    const auto *pfs_17 = buffer.data(pfs + 17);
    const auto *pfs_18 = buffer.data(pfs + 18);
    const auto *pfs_19 = buffer.data(pfs + 19);
    const auto *pfs_20 = buffer.data(pfs + 20);
    const auto *pfs_21 = buffer.data(pfs + 21);
    const auto *pfs_22 = buffer.data(pfs + 22);
    const auto *pfs_23 = buffer.data(pfs + 23);
    const auto *pfs_24 = buffer.data(pfs + 24);
    const auto *pfs_25 = buffer.data(pfs + 25);
    const auto *pfs_26 = buffer.data(pfs + 26);
    const auto *pfs_27 = buffer.data(pfs + 27);
    const auto *pfs_28 = buffer.data(pfs + 28);
    const auto *pfs_29 = buffer.data(pfs + 29);

    const auto *pfp1_0 = buffer.data(pfp1 + 0);
    const auto *pfp1_6 = buffer.data(pfp1 + 6);
    const auto *pfp1_9 = buffer.data(pfp1 + 9);
    const auto *pfp1_15 = buffer.data(pfp1 + 15);
    const auto *pfp1_34 = buffer.data(pfp1 + 34);
    const auto *pfp1_49 = buffer.data(pfp1 + 49);
    const auto *pfp1_52 = buffer.data(pfp1 + 52);
    const auto *pfp1_55 = buffer.data(pfp1 + 55);
    const auto *pfp1_68 = buffer.data(pfp1 + 68);
    const auto *pfp1_83 = buffer.data(pfp1 + 83);
    const auto *pfp1_86 = buffer.data(pfp1 + 86);
    const auto *pfp1_89 = buffer.data(pfp1 + 89);

    const auto *pgs_0 = buffer.data(pgs + 0);
    const auto *pgs_1 = buffer.data(pgs + 1);
    const auto *pgs_2 = buffer.data(pgs + 2);
    const auto *pgs_3 = buffer.data(pgs + 3);
    const auto *pgs_4 = buffer.data(pgs + 4);
    const auto *pgs_5 = buffer.data(pgs + 5);
    const auto *pgs_6 = buffer.data(pgs + 6);
    const auto *pgs_7 = buffer.data(pgs + 7);
    const auto *pgs_8 = buffer.data(pgs + 8);
    const auto *pgs_9 = buffer.data(pgs + 9);
    const auto *pgs_10 = buffer.data(pgs + 10);
    const auto *pgs_11 = buffer.data(pgs + 11);
    const auto *pgs_12 = buffer.data(pgs + 12);
    const auto *pgs_13 = buffer.data(pgs + 13);
    const auto *pgs_14 = buffer.data(pgs + 14);
    const auto *pgs_15 = buffer.data(pgs + 15);
    const auto *pgs_16 = buffer.data(pgs + 16);
    const auto *pgs_17 = buffer.data(pgs + 17);
    const auto *pgs_18 = buffer.data(pgs + 18);
    const auto *pgs_19 = buffer.data(pgs + 19);
    const auto *pgs_20 = buffer.data(pgs + 20);
    const auto *pgs_21 = buffer.data(pgs + 21);
    const auto *pgs_22 = buffer.data(pgs + 22);
    const auto *pgs_23 = buffer.data(pgs + 23);
    const auto *pgs_24 = buffer.data(pgs + 24);
    const auto *pgs_25 = buffer.data(pgs + 25);
    const auto *pgs_26 = buffer.data(pgs + 26);
    const auto *pgs_27 = buffer.data(pgs + 27);
    const auto *pgs_28 = buffer.data(pgs + 28);
    const auto *pgs_29 = buffer.data(pgs + 29);
    const auto *pgs_30 = buffer.data(pgs + 30);
    const auto *pgs_31 = buffer.data(pgs + 31);
    const auto *pgs_32 = buffer.data(pgs + 32);
    const auto *pgs_33 = buffer.data(pgs + 33);
    const auto *pgs_34 = buffer.data(pgs + 34);
    const auto *pgs_35 = buffer.data(pgs + 35);
    const auto *pgs_36 = buffer.data(pgs + 36);
    const auto *pgs_37 = buffer.data(pgs + 37);
    const auto *pgs_38 = buffer.data(pgs + 38);
    const auto *pgs_39 = buffer.data(pgs + 39);
    const auto *pgs_40 = buffer.data(pgs + 40);
    const auto *pgs_41 = buffer.data(pgs + 41);
    const auto *pgs_42 = buffer.data(pgs + 42);
    const auto *pgs_43 = buffer.data(pgs + 43);
    const auto *pgs_44 = buffer.data(pgs + 44);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_y, pc_x, pc_y, pc_z, sgs_0, pfp0_0, \
                         pfs_0, pfp1_0, pgs_0, pgs_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sgs_0[k]
                 + f_1 * pfs_0[k]
                 + f_2 * pc_x[k] * pgs_0[k];

        t_1[k] = f_2 * pc_y[k] * pgs_0[k];

        t_2[k] = f_2 * pc_z[k] * pgs_0[k];

        t_3[k] = pb_y[k] * pfp0_0[k]
                 - f_3 * pc_y[k] * pfp1_0[k];

        t_4[k] = f_0 * pfs_0[k]
                 + f_2 * pc_y[k] * pgs_1[k];

        t_5[k] = f_2 * pc_z[k] * pgs_1[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pb_z, pc_x, pc_y, pc_z, sgs_3, pfp0_0, pfs_0, \
                         pfs_3, pfp1_0, pgs_2, pgs_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pb_z[k] * pfp0_0[k]
                 - f_3 * pc_z[k] * pfp1_0[k];

        t_7[k] = f_2 * pc_y[k] * pgs_2[k];

        t_8[k] = f_0 * pfs_0[k]
                 + f_2 * pc_z[k] * pgs_2[k];

        t_9[k] = f_0 * sgs_3[k]
                 + f_4 * pfs_3[k]
                 + f_2 * pc_x[k] * pgs_3[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pb_y, pc_y, pc_z, pfp0_6, pfs_1, pfs_2, \
                         pfp1_6, pgs_3, pgs_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_4 * pfs_1[k]
                  + f_2 * pc_y[k] * pgs_3[k];

        t_11[k] = f_2 * pc_z[k] * pgs_3[k];

        t_12[k] = pb_y[k] * pfp0_6[k]
                  - f_3 * pc_y[k] * pfp1_6[k];

        t_13[k] = f_0 * pfs_2[k]
                  + f_2 * pc_y[k] * pgs_4[k];

        t_14[k] = f_0 * pfs_1[k]
                  + f_2 * pc_z[k] * pgs_4[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pc_x, pc_y, pc_z, sgs_5, sgs_6, pfs_2, \
                         pfs_3, pfs_5, pfs_6, pgs_5, pgs_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_0 * sgs_5[k]
                  + f_4 * pfs_5[k]
                  + f_2 * pc_x[k] * pgs_5[k];

        t_16[k] = f_2 * pc_y[k] * pgs_5[k];

        t_17[k] = f_4 * pfs_2[k]
                  + f_2 * pc_z[k] * pgs_5[k];

        t_18[k] = f_0 * sgs_6[k]
                  + f_0 * pfs_6[k]
                  + f_2 * pc_x[k] * pgs_6[k];

        t_19[k] = f_5 * pfs_3[k]
                  + f_2 * pc_y[k] * pgs_6[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pb_z, pc_y, pc_z, pfp0_9, pfs_3, pfs_4, \
                         pfp1_9, pgs_6, pgs_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_2 * pc_z[k] * pgs_6[k];

        t_21[k] = pb_z[k] * pfp0_9[k]
                  - f_3 * pc_z[k] * pfp1_9[k];

        t_22[k] = f_4 * pfs_4[k]
                  + f_2 * pc_y[k] * pgs_7[k];

        t_23[k] = f_0 * pfs_3[k]
                  + f_2 * pc_z[k] * pgs_7[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pb_y, pc_x, pc_y, pc_z, sgs_9, pfp0_15, \
                         pfs_4, pfs_5, pfs_9, pfp1_15, pgs_8, pgs_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pb_y[k] * pfp0_15[k]
                  - f_3 * pc_y[k] * pfp1_15[k];

        t_25[k] = f_0 * pfs_5[k]
                  + f_2 * pc_y[k] * pgs_8[k];

        t_26[k] = f_4 * pfs_4[k]
                  + f_2 * pc_z[k] * pgs_8[k];

        t_27[k] = f_0 * sgs_9[k]
                  + f_0 * pfs_9[k]
                  + f_2 * pc_x[k] * pgs_9[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, t_32, t_33, pc_x, pc_y, pc_z, sgs_10, sgs_11, \
                         pfs_5, pfs_6, pgs_9, pgs_10, pgs_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_2 * pc_y[k] * pgs_9[k];

        t_29[k] = f_5 * pfs_5[k]
                  + f_2 * pc_z[k] * pgs_9[k];

        t_30[k] = f_0 * sgs_10[k]
                  + f_2 * pc_x[k] * pgs_10[k];

        t_31[k] = f_1 * pfs_6[k]
                  + f_2 * pc_y[k] * pgs_10[k];

        t_32[k] = f_2 * pc_z[k] * pgs_10[k];

        t_33[k] = f_0 * sgs_11[k]
                  + f_2 * pc_x[k] * pgs_11[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, pc_x, pc_y, pc_z, sgs_12, pfs_6, pfs_7, \
                         pfs_8, pgs_11, pgs_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_5 * pfs_7[k]
                  + f_2 * pc_y[k] * pgs_11[k];

        t_35[k] = f_0 * pfs_6[k]
                  + f_2 * pc_z[k] * pgs_11[k];

        t_36[k] = f_0 * sgs_12[k]
                  + f_2 * pc_x[k] * pgs_12[k];

        t_37[k] = f_4 * pfs_8[k]
                  + f_2 * pc_y[k] * pgs_12[k];

        t_38[k] = f_4 * pfs_7[k]
                  + f_2 * pc_z[k] * pgs_12[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, t_44, pc_x, pc_y, pc_z, sgs_13, sgs_14, \
                         pfs_8, pfs_9, pgs_13, pgs_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_0 * sgs_13[k]
                  + f_2 * pc_x[k] * pgs_13[k];

        t_40[k] = f_0 * pfs_9[k]
                  + f_2 * pc_y[k] * pgs_13[k];

        t_41[k] = f_5 * pfs_8[k]
                  + f_2 * pc_z[k] * pgs_13[k];

        t_42[k] = f_0 * sgs_14[k]
                  + f_2 * pc_x[k] * pgs_14[k];

        t_43[k] = f_2 * pc_y[k] * pgs_14[k];

        t_44[k] = f_1 * pfs_9[k]
                  + f_2 * pc_z[k] * pgs_14[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, t_50, pc_x, pc_y, pc_z, sgs_0, sgs_1, \
                         pfs_10, pfs_11, pgs_15, pgs_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_1 * pfs_10[k]
                  + f_2 * pc_x[k] * pgs_15[k];

        t_46[k] = f_0 * sgs_0[k]
                  + f_2 * pc_y[k] * pgs_15[k];

        t_47[k] = f_2 * pc_z[k] * pgs_15[k];

        t_48[k] = f_5 * pfs_11[k]
                  + f_2 * pc_x[k] * pgs_16[k];

        t_49[k] = f_0 * sgs_1[k]
                  + f_0 * pfs_10[k]
                  + f_2 * pc_y[k] * pgs_16[k];

        t_50[k] = f_2 * pc_z[k] * pgs_16[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, pc_x, pc_y, pc_z, sgs_2, sgs_3, pfs_10, \
                         pfs_11, pfs_12, pfs_13, pgs_17, pgs_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_5 * pfs_12[k]
                  + f_2 * pc_x[k] * pgs_17[k];

        t_52[k] = f_0 * sgs_2[k]
                  + f_2 * pc_y[k] * pgs_17[k];

        t_53[k] = f_0 * pfs_10[k]
                  + f_2 * pc_z[k] * pgs_17[k];

        t_54[k] = f_4 * pfs_13[k]
                  + f_2 * pc_x[k] * pgs_18[k];

        t_55[k] = f_0 * sgs_3[k]
                  + f_4 * pfs_11[k]
                  + f_2 * pc_y[k] * pgs_18[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, pb_z, pc_x, pc_z, pfp0_34, pfs_11, \
                         pfs_14, pfs_15, pfp1_34, pgs_18, pgs_19, \
                         pgs_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_2 * pc_z[k] * pgs_18[k];

        t_57[k] = f_4 * pfs_14[k]
                  + f_2 * pc_x[k] * pgs_19[k];

        t_58[k] = pb_z[k] * pfp0_34[k]
                  - f_3 * pc_z[k] * pfp1_34[k];

        t_59[k] = f_0 * pfs_11[k]
                  + f_2 * pc_z[k] * pgs_19[k];

        t_60[k] = f_4 * pfs_15[k]
                  + f_2 * pc_x[k] * pgs_20[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, t_65, pb_x, pc_x, pc_y, pc_z, sgs_5, pfp0_49, \
                         pfs_12, pfs_16, pfp1_49, pgs_20, pgs_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_0 * sgs_5[k]
                  + f_2 * pc_y[k] * pgs_20[k];

        t_62[k] = f_4 * pfs_12[k]
                  + f_2 * pc_z[k] * pgs_20[k];

        t_63[k] = f_0 * pfs_16[k]
                  + f_2 * pc_x[k] * pgs_21[k];

        t_64[k] = pb_x[k] * pfp0_49[k]
                  - f_3 * pc_x[k] * pfp1_49[k];

        t_65[k] = f_2 * pc_z[k] * pgs_21[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pb_x, pc_x, pc_z, pfp0_52, pfs_13, pfs_17, \
                         pfs_18, pfp1_52, pgs_22, pgs_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_0 * pfs_17[k]
                  + f_2 * pc_x[k] * pgs_22[k];

        t_67[k] = pb_x[k] * pfp0_52[k]
                  - f_3 * pc_x[k] * pfp1_52[k];

        t_68[k] = f_0 * pfs_13[k]
                  + f_2 * pc_z[k] * pgs_22[k];

        t_69[k] = f_0 * pfs_18[k]
                  + f_2 * pc_x[k] * pgs_23[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pb_x, pc_x, pc_y, pc_z, sgs_9, pfp0_55, \
                         pfs_14, pfs_19, pfp1_55, pgs_23, pgs_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pb_x[k] * pfp0_55[k]
                  - f_3 * pc_x[k] * pfp1_55[k];

        t_71[k] = f_4 * pfs_14[k]
                  + f_2 * pc_z[k] * pgs_23[k];

        t_72[k] = f_0 * pfs_19[k]
                  + f_2 * pc_x[k] * pgs_24[k];

        t_73[k] = f_0 * sgs_9[k]
                  + f_2 * pc_y[k] * pgs_24[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, t_78, pc_x, pc_y, pc_z, sgs_10, pfs_15, \
                         pfs_16, pgs_24, pgs_25, pgs_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_5 * pfs_15[k]
                  + f_2 * pc_z[k] * pgs_24[k];

        t_75[k] = f_2 * pc_x[k] * pgs_25[k];

        t_76[k] = f_0 * sgs_10[k]
                  + f_1 * pfs_16[k]
                  + f_2 * pc_y[k] * pgs_25[k];

        t_77[k] = f_2 * pc_z[k] * pgs_25[k];

        t_78[k] = f_2 * pc_x[k] * pgs_26[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pb_z, pc_x, pc_y, pc_z, sgs_12, pfp0_49, \
                         pfs_16, pfs_18, pfp1_49, pgs_26, pgs_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = pb_z[k] * pfp0_49[k]
                  - f_3 * pc_z[k] * pfp1_49[k];

        t_80[k] = f_0 * pfs_16[k]
                  + f_2 * pc_z[k] * pgs_26[k];

        t_81[k] = f_2 * pc_x[k] * pgs_27[k];

        t_82[k] = f_0 * sgs_12[k]
                  + f_4 * pfs_18[k]
                  + f_2 * pc_y[k] * pgs_27[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, t_87, pc_x, pc_y, pc_z, sgs_13, pfs_17, \
                         pfs_18, pfs_19, pgs_27, pgs_28, pgs_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_4 * pfs_17[k]
                  + f_2 * pc_z[k] * pgs_27[k];

        t_84[k] = f_2 * pc_x[k] * pgs_28[k];

        t_85[k] = f_0 * sgs_13[k]
                  + f_0 * pfs_19[k]
                  + f_2 * pc_y[k] * pgs_28[k];

        t_86[k] = f_5 * pfs_18[k]
                  + f_2 * pc_z[k] * pgs_28[k];

        t_87[k] = f_2 * pc_x[k] * pgs_29[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, t_92, pc_x, pc_y, pc_z, sgs_0, sgs_14, \
                         pfs_19, pfs_20, pgs_29, pgs_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_0 * sgs_14[k]
                  + f_2 * pc_y[k] * pgs_29[k];

        t_89[k] = f_1 * pfs_19[k]
                  + f_2 * pc_z[k] * pgs_29[k];

        t_90[k] = f_1 * pfs_20[k]
                  + f_2 * pc_x[k] * pgs_30[k];

        t_91[k] = f_2 * pc_y[k] * pgs_30[k];

        t_92[k] = f_0 * sgs_0[k]
                  + f_2 * pc_z[k] * pgs_30[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, t_97, t_98, pc_x, pc_y, pc_z, sgs_1, sgs_2, \
                         pfs_20, pfs_21, pfs_22, pgs_31, pgs_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_5 * pfs_21[k]
                  + f_2 * pc_x[k] * pgs_31[k];

        t_94[k] = f_0 * pfs_20[k]
                  + f_2 * pc_y[k] * pgs_31[k];

        t_95[k] = f_0 * sgs_1[k]
                  + f_2 * pc_z[k] * pgs_31[k];

        t_96[k] = f_5 * pfs_22[k]
                  + f_2 * pc_x[k] * pgs_32[k];

        t_97[k] = f_2 * pc_y[k] * pgs_32[k];

        t_98[k] = f_0 * sgs_2[k]
                  + f_0 * pfs_20[k]
                  + f_2 * pc_z[k] * pgs_32[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, t_103, pc_x, pc_y, pc_z, sgs_3, pfs_21, \
                         pfs_22, pfs_23, pfs_24, pgs_33, pgs_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_4 * pfs_23[k]
                  + f_2 * pc_x[k] * pgs_33[k];

        t_100[k] = f_4 * pfs_21[k]
                   + f_2 * pc_y[k] * pgs_33[k];

        t_101[k] = f_0 * sgs_3[k]
                   + f_2 * pc_z[k] * pgs_33[k];

        t_102[k] = f_4 * pfs_24[k]
                   + f_2 * pc_x[k] * pgs_34[k];

        t_103[k] = f_0 * pfs_22[k]
                   + f_2 * pc_y[k] * pgs_34[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pb_y, pc_x, pc_y, pc_z, sgs_5, pfp0_68, \
                         pfs_22, pfs_25, pfp1_68, pgs_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pb_y[k] * pfp0_68[k]
                   - f_3 * pc_y[k] * pfp1_68[k];

        t_105[k] = f_4 * pfs_25[k]
                   + f_2 * pc_x[k] * pgs_35[k];

        t_106[k] = f_2 * pc_y[k] * pgs_35[k];

        t_107[k] = f_0 * sgs_5[k]
                   + f_4 * pfs_22[k]
                   + f_2 * pc_z[k] * pgs_35[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, t_111, t_112, pc_x, pc_y, pc_z, sgs_6, pfs_23, \
                         pfs_24, pfs_26, pfs_27, pgs_36, pgs_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_0 * pfs_26[k]
                   + f_2 * pc_x[k] * pgs_36[k];

        t_109[k] = f_5 * pfs_23[k]
                   + f_2 * pc_y[k] * pgs_36[k];

        t_110[k] = f_0 * sgs_6[k]
                   + f_2 * pc_z[k] * pgs_36[k];

        t_111[k] = f_0 * pfs_27[k]
                   + f_2 * pc_x[k] * pgs_37[k];

        t_112[k] = f_4 * pfs_24[k]
                   + f_2 * pc_y[k] * pgs_37[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pb_x, pc_x, pc_y, pfp0_83, pfp0_86, \
                         pfs_25, pfs_28, pfp1_83, pfp1_86, pgs_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = pb_x[k] * pfp0_83[k]
                   - f_3 * pc_x[k] * pfp1_83[k];

        t_114[k] = f_0 * pfs_28[k]
                   + f_2 * pc_x[k] * pgs_38[k];

        t_115[k] = f_0 * pfs_25[k]
                   + f_2 * pc_y[k] * pgs_38[k];

        t_116[k] = pb_x[k] * pfp0_86[k]
                   - f_3 * pc_x[k] * pfp1_86[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, t_121, pb_x, pc_x, pc_y, pfp0_89, pfs_26, \
                         pfs_29, pfp1_89, pgs_39, pgs_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_0 * pfs_29[k]
                   + f_2 * pc_x[k] * pgs_39[k];

        t_118[k] = f_2 * pc_y[k] * pgs_39[k];

        t_119[k] = pb_x[k] * pfp0_89[k]
                   - f_3 * pc_x[k] * pfp1_89[k];

        t_120[k] = f_2 * pc_x[k] * pgs_40[k];

        t_121[k] = f_1 * pfs_26[k]
                   + f_2 * pc_y[k] * pgs_40[k];
    }

#pragma omp simd aligned(t_122, t_123, t_124, t_125, t_126, pc_x, pc_y, pc_z, sgs_10, sgs_11, \
                         pfs_26, pfs_27, pgs_40, pgs_41, pgs_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_122[k] = f_0 * sgs_10[k]
                   + f_2 * pc_z[k] * pgs_40[k];

        t_123[k] = f_2 * pc_x[k] * pgs_41[k];

        t_124[k] = f_5 * pfs_27[k]
                   + f_2 * pc_y[k] * pgs_41[k];

        t_125[k] = f_0 * sgs_11[k]
                   + f_0 * pfs_26[k]
                   + f_2 * pc_z[k] * pgs_41[k];

        t_126[k] = f_2 * pc_x[k] * pgs_42[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, pc_x, pc_y, pc_z, sgs_12, pfs_27, pfs_28, \
                         pfs_29, pgs_42, pgs_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_4 * pfs_28[k]
                   + f_2 * pc_y[k] * pgs_42[k];

        t_128[k] = f_0 * sgs_12[k]
                   + f_4 * pfs_27[k]
                   + f_2 * pc_z[k] * pgs_42[k];

        t_129[k] = f_2 * pc_x[k] * pgs_43[k];

        t_130[k] = f_0 * pfs_29[k]
                   + f_2 * pc_y[k] * pgs_43[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, pb_y, pc_x, pc_y, pc_z, sgs_14, pfp0_89, \
                         pfs_29, pfp1_89, pgs_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = pb_y[k] * pfp0_89[k]
                   - f_3 * pc_y[k] * pfp1_89[k];

        t_132[k] = f_2 * pc_x[k] * pgs_44[k];

        t_133[k] = f_2 * pc_y[k] * pgs_44[k];

        t_134[k] = f_0 * sgs_14[k]
                   + f_1 * pfs_29[k]
                   + f_2 * pc_z[k] * pgs_44[k];
    }
}

}  // namespace simdt3ceri
