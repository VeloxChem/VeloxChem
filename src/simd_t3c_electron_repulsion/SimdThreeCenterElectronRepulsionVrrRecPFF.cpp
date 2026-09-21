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


#include "SimdThreeCenterElectronRepulsionVrrRecPFF.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_pff_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sff0, const size_t sfd,
                                                          const size_t sff1, const size_t ppf0,
                                                          const size_t ppf1, const size_t pdf0,
                                                          const size_t pdd, const size_t pdf1,
                                                          const size_t pfp0, const size_t pfp1,
                                                          const size_t pfd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.5 / q;
    const auto f_2 = 1.0 / gamma;
    const auto f_3 = p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = gamma / q;
    const auto f_6 = 1.0 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_9 = 0.5 / p;
    const auto f_10 = 0.5 * gamma / (p * q);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sff0_0 = buffer.data(sff0 + 0);
    const auto *sff0_1 = buffer.data(sff0 + 1);
    const auto *sff0_6 = buffer.data(sff0 + 6);
    const auto *sff0_9 = buffer.data(sff0 + 9);
    const auto *sff0_20 = buffer.data(sff0 + 20);
    const auto *sff0_60 = buffer.data(sff0 + 60);
    const auto *sff0_66 = buffer.data(sff0 + 66);
    const auto *sff0_69 = buffer.data(sff0 + 69);
    const auto *sff0_76 = buffer.data(sff0 + 76);
    const auto *sff0_79 = buffer.data(sff0 + 79);
    const auto *sff0_86 = buffer.data(sff0 + 86);
    const auto *sff0_89 = buffer.data(sff0 + 89);
    const auto *sff0_90 = buffer.data(sff0 + 90);
    const auto *sff0_96 = buffer.data(sff0 + 96);
    const auto *sff0_99 = buffer.data(sff0 + 99);

    const auto *sfd_0 = buffer.data(sfd + 0);
    const auto *sfd_3 = buffer.data(sfd + 3);
    const auto *sfd_5 = buffer.data(sfd + 5);
    const auto *sfd_9 = buffer.data(sfd + 9);
    const auto *sfd_11 = buffer.data(sfd + 11);
    const auto *sfd_17 = buffer.data(sfd + 17);
    const auto *sfd_18 = buffer.data(sfd + 18);
    const auto *sfd_21 = buffer.data(sfd + 21);
    const auto *sfd_23 = buffer.data(sfd + 23);
    const auto *sfd_30 = buffer.data(sfd + 30);
    const auto *sfd_33 = buffer.data(sfd + 33);
    const auto *sfd_35 = buffer.data(sfd + 35);
    const auto *sfd_36 = buffer.data(sfd + 36);
    const auto *sfd_39 = buffer.data(sfd + 39);
    const auto *sfd_41 = buffer.data(sfd + 41);
    const auto *sfd_45 = buffer.data(sfd + 45);
    const auto *sfd_47 = buffer.data(sfd + 47);
    const auto *sfd_51 = buffer.data(sfd + 51);
    const auto *sfd_53 = buffer.data(sfd + 53);
    const auto *sfd_54 = buffer.data(sfd + 54);
    const auto *sfd_57 = buffer.data(sfd + 57);
    const auto *sfd_59 = buffer.data(sfd + 59);

    const auto *sff1_0 = buffer.data(sff1 + 0);
    const auto *sff1_1 = buffer.data(sff1 + 1);
    const auto *sff1_6 = buffer.data(sff1 + 6);
    const auto *sff1_9 = buffer.data(sff1 + 9);
    const auto *sff1_20 = buffer.data(sff1 + 20);
    const auto *sff1_60 = buffer.data(sff1 + 60);
    const auto *sff1_66 = buffer.data(sff1 + 66);
    const auto *sff1_69 = buffer.data(sff1 + 69);
    const auto *sff1_76 = buffer.data(sff1 + 76);
    const auto *sff1_79 = buffer.data(sff1 + 79);
    const auto *sff1_86 = buffer.data(sff1 + 86);
    const auto *sff1_89 = buffer.data(sff1 + 89);
    const auto *sff1_90 = buffer.data(sff1 + 90);
    const auto *sff1_96 = buffer.data(sff1 + 96);
    const auto *sff1_99 = buffer.data(sff1 + 99);

    const auto *ppf0_46 = buffer.data(ppf0 + 46);

    const auto *ppf1_46 = buffer.data(ppf1 + 46);

    const auto *pdf0_0 = buffer.data(pdf0 + 0);
    const auto *pdf0_3 = buffer.data(pdf0 + 3);
    const auto *pdf0_5 = buffer.data(pdf0 + 5);
    const auto *pdf0_6 = buffer.data(pdf0 + 6);
    const auto *pdf0_9 = buffer.data(pdf0 + 9);
    const auto *pdf0_13 = buffer.data(pdf0 + 13);
    const auto *pdf0_20 = buffer.data(pdf0 + 20);
    const auto *pdf0_25 = buffer.data(pdf0 + 25);
    const auto *pdf0_30 = buffer.data(pdf0 + 30);
    const auto *pdf0_50 = buffer.data(pdf0 + 50);
    const auto *pdf0_61 = buffer.data(pdf0 + 61);
    const auto *pdf0_66 = buffer.data(pdf0 + 66);
    const auto *pdf0_76 = buffer.data(pdf0 + 76);

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
    const auto *pdd_43 = buffer.data(pdd + 43);
    const auto *pdd_45 = buffer.data(pdd + 45);
    const auto *pdd_46 = buffer.data(pdd + 46);
    const auto *pdd_47 = buffer.data(pdd + 47);
    const auto *pdd_51 = buffer.data(pdd + 51);
    const auto *pdd_52 = buffer.data(pdd + 52);
    const auto *pdd_53 = buffer.data(pdd + 53);

    const auto *pdf1_0 = buffer.data(pdf1 + 0);
    const auto *pdf1_3 = buffer.data(pdf1 + 3);
    const auto *pdf1_5 = buffer.data(pdf1 + 5);
    const auto *pdf1_6 = buffer.data(pdf1 + 6);
    const auto *pdf1_9 = buffer.data(pdf1 + 9);
    const auto *pdf1_13 = buffer.data(pdf1 + 13);
    const auto *pdf1_20 = buffer.data(pdf1 + 20);
    const auto *pdf1_25 = buffer.data(pdf1 + 25);
    const auto *pdf1_30 = buffer.data(pdf1 + 30);
    const auto *pdf1_50 = buffer.data(pdf1 + 50);
    const auto *pdf1_61 = buffer.data(pdf1 + 61);
    const auto *pdf1_66 = buffer.data(pdf1 + 66);
    const auto *pdf1_76 = buffer.data(pdf1 + 76);

    const auto *pfp0_0 = buffer.data(pfp0 + 0);
    const auto *pfp0_1 = buffer.data(pfp0 + 1);
    const auto *pfp0_2 = buffer.data(pfp0 + 2);
    const auto *pfp0_4 = buffer.data(pfp0 + 4);
    const auto *pfp0_8 = buffer.data(pfp0 + 8);
    const auto *pfp0_9 = buffer.data(pfp0 + 9);
    const auto *pfp0_10 = buffer.data(pfp0 + 10);
    const auto *pfp0_11 = buffer.data(pfp0 + 11);
    const auto *pfp0_13 = buffer.data(pfp0 + 13);
    const auto *pfp0_14 = buffer.data(pfp0 + 14);
    const auto *pfp0_15 = buffer.data(pfp0 + 15);
    const auto *pfp0_16 = buffer.data(pfp0 + 16);
    const auto *pfp0_17 = buffer.data(pfp0 + 17);
    const auto *pfp0_33 = buffer.data(pfp0 + 33);
    const auto *pfp0_34 = buffer.data(pfp0 + 34);
    const auto *pfp0_35 = buffer.data(pfp0 + 35);

    const auto *pfp1_0 = buffer.data(pfp1 + 0);
    const auto *pfp1_1 = buffer.data(pfp1 + 1);
    const auto *pfp1_2 = buffer.data(pfp1 + 2);
    const auto *pfp1_4 = buffer.data(pfp1 + 4);
    const auto *pfp1_8 = buffer.data(pfp1 + 8);
    const auto *pfp1_9 = buffer.data(pfp1 + 9);
    const auto *pfp1_10 = buffer.data(pfp1 + 10);
    const auto *pfp1_11 = buffer.data(pfp1 + 11);
    const auto *pfp1_13 = buffer.data(pfp1 + 13);
    const auto *pfp1_14 = buffer.data(pfp1 + 14);
    const auto *pfp1_15 = buffer.data(pfp1 + 15);
    const auto *pfp1_16 = buffer.data(pfp1 + 16);
    const auto *pfp1_17 = buffer.data(pfp1 + 17);
    const auto *pfp1_33 = buffer.data(pfp1 + 33);
    const auto *pfp1_34 = buffer.data(pfp1 + 34);
    const auto *pfp1_35 = buffer.data(pfp1 + 35);

    const auto *pfd_0 = buffer.data(pfd + 0);
    const auto *pfd_2 = buffer.data(pfd + 2);
    const auto *pfd_3 = buffer.data(pfd + 3);
    const auto *pfd_5 = buffer.data(pfd + 5);
    const auto *pfd_6 = buffer.data(pfd + 6);
    const auto *pfd_8 = buffer.data(pfd + 8);
    const auto *pfd_9 = buffer.data(pfd + 9);
    const auto *pfd_11 = buffer.data(pfd + 11);
    const auto *pfd_12 = buffer.data(pfd + 12);
    const auto *pfd_14 = buffer.data(pfd + 14);
    const auto *pfd_15 = buffer.data(pfd + 15);
    const auto *pfd_17 = buffer.data(pfd + 17);
    const auto *pfd_18 = buffer.data(pfd + 18);
    const auto *pfd_20 = buffer.data(pfd + 20);
    const auto *pfd_21 = buffer.data(pfd + 21);
    const auto *pfd_23 = buffer.data(pfd + 23);
    const auto *pfd_24 = buffer.data(pfd + 24);
    const auto *pfd_26 = buffer.data(pfd + 26);
    const auto *pfd_27 = buffer.data(pfd + 27);
    const auto *pfd_29 = buffer.data(pfd + 29);
    const auto *pfd_30 = buffer.data(pfd + 30);
    const auto *pfd_32 = buffer.data(pfd + 32);
    const auto *pfd_33 = buffer.data(pfd + 33);
    const auto *pfd_35 = buffer.data(pfd + 35);
    const auto *pfd_36 = buffer.data(pfd + 36);
    const auto *pfd_38 = buffer.data(pfd + 38);
    const auto *pfd_39 = buffer.data(pfd + 39);
    const auto *pfd_41 = buffer.data(pfd + 41);
    const auto *pfd_42 = buffer.data(pfd + 42);
    const auto *pfd_44 = buffer.data(pfd + 44);
    const auto *pfd_45 = buffer.data(pfd + 45);
    const auto *pfd_47 = buffer.data(pfd + 47);
    const auto *pfd_48 = buffer.data(pfd + 48);
    const auto *pfd_50 = buffer.data(pfd + 50);
    const auto *pfd_51 = buffer.data(pfd + 51);
    const auto *pfd_53 = buffer.data(pfd + 53);
    const auto *pfd_54 = buffer.data(pfd + 54);
    const auto *pfd_56 = buffer.data(pfd + 56);
    const auto *pfd_57 = buffer.data(pfd + 57);
    const auto *pfd_59 = buffer.data(pfd + 59);
    const auto *pfd_60 = buffer.data(pfd + 60);
    const auto *pfd_63 = buffer.data(pfd + 63);
    const auto *pfd_64 = buffer.data(pfd + 64);
    const auto *pfd_65 = buffer.data(pfd + 65);
    const auto *pfd_66 = buffer.data(pfd + 66);
    const auto *pfd_67 = buffer.data(pfd + 67);
    const auto *pfd_69 = buffer.data(pfd + 69);
    const auto *pfd_70 = buffer.data(pfd + 70);
    const auto *pfd_71 = buffer.data(pfd + 71);
    const auto *pfd_72 = buffer.data(pfd + 72);
    const auto *pfd_75 = buffer.data(pfd + 75);
    const auto *pfd_76 = buffer.data(pfd + 76);
    const auto *pfd_77 = buffer.data(pfd + 77);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, sfd_0, sfd_3, pdd_0, pdd_3, \
                         pfp0_0, pfp1_0, pfd_0, pfd_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sfd_0[k]
                 + f_1 * pdd_0[k]
                 + f_2 * pfp0_0[k]
                 - f_3 * pfp1_0[k]
                 + f_4 * pc_x[k] * pfd_0[k];

        t_1[k] = f_4 * pc_y[k] * pfd_0[k];

        t_2[k] = f_4 * pc_z[k] * pfd_0[k];

        t_3[k] = f_0 * sfd_3[k]
                 + f_1 * pdd_3[k]
                 + f_4 * pc_x[k] * pfd_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, t_8, pc_x, pc_y, pc_z, sfd_5, pdd_5, pfp0_1, \
                         pfp1_1, pfd_2, pfd_3, pfd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_4 * pc_y[k] * pfd_2[k];

        t_5[k] = f_0 * sfd_5[k]
                 + f_1 * pdd_5[k]
                 + f_4 * pc_x[k] * pfd_5[k];

        t_6[k] = f_2 * pfp0_1[k]
                 - f_3 * pfp1_1[k]
                 + f_4 * pc_y[k] * pfd_3[k];

        t_7[k] = f_4 * pc_z[k] * pfd_3[k];

        t_8[k] = f_4 * pc_y[k] * pfd_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pb_y, pc_y, pc_z, pdf0_0, pdd_0, pdf1_0, \
                         pfp0_2, pfp1_2, pfd_5, pfd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_2 * pfp0_2[k]
                 - f_3 * pfp1_2[k]
                 + f_4 * pc_z[k] * pfd_5[k];

        t_10[k] = pb_y[k] * pdf0_0[k]
                  - f_5 * pc_y[k] * pdf1_0[k];

        t_11[k] = f_0 * pdd_0[k]
                  + f_4 * pc_y[k] * pfd_6[k];

        t_12[k] = f_4 * pc_z[k] * pfd_6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pb_y, pc_x, pc_y, sfd_9, pdf0_5, pdd_2, pdd_9, \
                         pdf1_5, pfd_8, pfd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_0 * sfd_9[k]
                  + f_6 * pdd_9[k]
                  + f_4 * pc_x[k] * pfd_9[k];

        t_14[k] = f_0 * pdd_2[k]
                  + f_4 * pc_y[k] * pfd_8[k];

        t_15[k] = pb_y[k] * pdf0_5[k]
                  - f_5 * pc_y[k] * pdf1_5[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pb_y, pc_y, pc_z, pdf0_9, pdd_3, pdd_5, \
                         pdf1_9, pfp0_4, pfp1_4, pfd_9, pfd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_0 * pdd_3[k]
                  + f_2 * pfp0_4[k]
                  - f_3 * pfp1_4[k]
                  + f_4 * pc_y[k] * pfd_9[k];

        t_17[k] = f_4 * pc_z[k] * pfd_9[k];

        t_18[k] = f_0 * pdd_5[k]
                  + f_4 * pc_y[k] * pfd_11[k];

        t_19[k] = pb_y[k] * pdf0_9[k]
                  - f_5 * pc_y[k] * pdf1_9[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pb_z, pc_y, pc_z, pdf0_0, pdf0_3, \
                         pdd_0, pdf1_0, pdf1_3, pfd_12, pfd_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pb_z[k] * pdf0_0[k]
                  - f_5 * pc_z[k] * pdf1_0[k];

        t_21[k] = f_4 * pc_y[k] * pfd_12[k];

        t_22[k] = f_0 * pdd_0[k]
                  + f_4 * pc_z[k] * pfd_12[k];

        t_23[k] = pb_z[k] * pdf0_3[k]
                  - f_5 * pc_z[k] * pdf1_3[k];

        t_24[k] = f_4 * pc_y[k] * pfd_14[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pb_z, pc_x, pc_y, pc_z, sfd_17, pdf0_6, \
                         pdd_3, pdd_17, pdf1_6, pfd_15, pfd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_0 * sfd_17[k]
                  + f_6 * pdd_17[k]
                  + f_4 * pc_x[k] * pfd_17[k];

        t_26[k] = pb_z[k] * pdf0_6[k]
                  - f_5 * pc_z[k] * pdf1_6[k];

        t_27[k] = f_0 * pdd_3[k]
                  + f_4 * pc_z[k] * pfd_15[k];

        t_28[k] = f_4 * pc_y[k] * pfd_17[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, pc_x, pc_y, pc_z, sfd_18, pdd_5, pdd_6, pdd_18, \
                         pfp0_8, pfp0_9, pfp1_8, pfp1_9, pfd_17, \
                         pfd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * pdd_5[k]
                  + f_2 * pfp0_8[k]
                  - f_3 * pfp1_8[k]
                  + f_4 * pc_z[k] * pfd_17[k];

        t_30[k] = f_0 * sfd_18[k]
                  + f_0 * pdd_18[k]
                  + f_2 * pfp0_9[k]
                  - f_3 * pfp1_9[k]
                  + f_4 * pc_x[k] * pfd_18[k];

        t_31[k] = f_6 * pdd_6[k]
                  + f_4 * pc_y[k] * pfd_18[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pc_x, pc_y, pc_z, sfd_21, sfd_23, pdd_8, \
                         pdd_21, pdd_23, pfd_18, pfd_20, pfd_21, \
                         pfd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_4 * pc_z[k] * pfd_18[k];

        t_33[k] = f_0 * sfd_21[k]
                  + f_0 * pdd_21[k]
                  + f_4 * pc_x[k] * pfd_21[k];

        t_34[k] = f_6 * pdd_8[k]
                  + f_4 * pc_y[k] * pfd_20[k];

        t_35[k] = f_0 * sfd_23[k]
                  + f_0 * pdd_23[k]
                  + f_4 * pc_x[k] * pfd_23[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pc_y, pc_z, pdd_9, pdd_11, pfp0_10, pfp0_11, \
                         pfp1_10, pfp1_11, pfd_21, pfd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_6 * pdd_9[k]
                  + f_2 * pfp0_10[k]
                  - f_3 * pfp1_10[k]
                  + f_4 * pc_y[k] * pfd_21[k];

        t_37[k] = f_4 * pc_z[k] * pfd_21[k];

        t_38[k] = f_6 * pdd_11[k]
                  + f_4 * pc_y[k] * pfd_23[k];

        t_39[k] = f_2 * pfp0_11[k]
                  - f_3 * pfp1_11[k]
                  + f_4 * pc_z[k] * pfd_23[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pb_y, pb_z, pc_y, pc_z, pdf0_13, pdf0_20, \
                         pdd_6, pdd_12, pdf1_13, pdf1_20, pfd_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_y[k] * pdf0_20[k]
                  - f_5 * pc_y[k] * pdf1_20[k];

        t_41[k] = f_0 * pdd_12[k]
                  + f_4 * pc_y[k] * pfd_24[k];

        t_42[k] = f_0 * pdd_6[k]
                  + f_4 * pc_z[k] * pfd_24[k];

        t_43[k] = pb_z[k] * pdf0_13[k]
                  - f_5 * pc_z[k] * pdf1_13[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pb_y, pc_y, pc_z, pdf0_25, pdd_9, pdd_14, \
                         pdd_15, pdf1_25, pfp0_13, pfp1_13, pfd_26, \
                         pfd_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_0 * pdd_14[k]
                  + f_4 * pc_y[k] * pfd_26[k];

        t_45[k] = pb_y[k] * pdf0_25[k]
                  - f_5 * pc_y[k] * pdf1_25[k];

        t_46[k] = f_0 * pdd_15[k]
                  + f_2 * pfp0_13[k]
                  - f_3 * pfp1_13[k]
                  + f_4 * pc_y[k] * pfd_27[k];

        t_47[k] = f_0 * pdd_9[k]
                  + f_4 * pc_z[k] * pfd_27[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, pc_x, pc_y, pc_z, sfd_30, pdd_11, pdd_17, pdd_30, \
                         pfp0_14, pfp0_15, pfp1_14, pfp1_15, pfd_29, \
                         pfd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_0 * pdd_17[k]
                  + f_4 * pc_y[k] * pfd_29[k];

        t_49[k] = f_0 * pdd_11[k]
                  + f_2 * pfp0_14[k]
                  - f_3 * pfp1_14[k]
                  + f_4 * pc_z[k] * pfd_29[k];

        t_50[k] = f_0 * sfd_30[k]
                  + f_0 * pdd_30[k]
                  + f_2 * pfp0_15[k]
                  - f_3 * pfp1_15[k]
                  + f_4 * pc_x[k] * pfd_30[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pc_x, pc_y, pc_z, sfd_33, pdd_12, pdd_33, \
                         pfd_30, pfd_32, pfd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_4 * pc_y[k] * pfd_30[k];

        t_52[k] = f_6 * pdd_12[k]
                  + f_4 * pc_z[k] * pfd_30[k];

        t_53[k] = f_0 * sfd_33[k]
                  + f_0 * pdd_33[k]
                  + f_4 * pc_x[k] * pfd_33[k];

        t_54[k] = f_4 * pc_y[k] * pfd_32[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pc_x, pc_y, pc_z, sfd_35, pdd_15, pdd_35, \
                         pfp0_16, pfp1_16, pfd_33, pfd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_0 * sfd_35[k]
                  + f_0 * pdd_35[k]
                  + f_4 * pc_x[k] * pfd_35[k];

        t_56[k] = f_2 * pfp0_16[k]
                  - f_3 * pfp1_16[k]
                  + f_4 * pc_y[k] * pfd_33[k];

        t_57[k] = f_6 * pdd_15[k]
                  + f_4 * pc_z[k] * pfd_33[k];

        t_58[k] = f_4 * pc_y[k] * pfd_35[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, pa_x, pc_x, pc_y, pc_z, sff0_60, sfd_36, sff1_60, \
                         pdd_17, pdd_18, pfp0_17, pfp1_17, pfd_35, \
                         pfd_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_6 * pdd_17[k]
                  + f_2 * pfp0_17[k]
                  - f_3 * pfp1_17[k]
                  + f_4 * pc_z[k] * pfd_35[k];

        t_60[k] = pa_x[k] * sff0_60[k]
                  + f_1 * sfd_36[k]
                  - f_5 * pc_x[k] * sff1_60[k];

        t_61[k] = f_1 * pdd_18[k]
                  + f_4 * pc_y[k] * pfd_36[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pc_x, pc_y, pc_z, sfd_39, sfd_41, pdd_20, \
                         pfd_36, pfd_38, pfd_39, pfd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_4 * pc_z[k] * pfd_36[k];

        t_63[k] = f_0 * sfd_39[k]
                  + f_4 * pc_x[k] * pfd_39[k];

        t_64[k] = f_1 * pdd_20[k]
                  + f_4 * pc_y[k] * pfd_38[k];

        t_65[k] = f_0 * sfd_41[k]
                  + f_4 * pc_x[k] * pfd_41[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pa_x, pc_x, pc_y, pc_z, sff0_66, sff0_69, \
                         sff1_66, sff1_69, pdd_23, pfd_39, pfd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pa_x[k] * sff0_66[k]
                  - f_5 * pc_x[k] * sff1_66[k];

        t_67[k] = f_4 * pc_z[k] * pfd_39[k];

        t_68[k] = f_1 * pdd_23[k]
                  + f_4 * pc_y[k] * pfd_41[k];

        t_69[k] = pa_x[k] * sff0_69[k]
                  - f_5 * pc_x[k] * sff1_69[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pb_z, pc_x, pc_y, pc_z, sfd_45, pdf0_30, \
                         pdd_18, pdd_24, pdf1_30, pfd_42, pfd_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pb_z[k] * pdf0_30[k]
                  - f_5 * pc_z[k] * pdf1_30[k];

        t_71[k] = f_6 * pdd_24[k]
                  + f_4 * pc_y[k] * pfd_42[k];

        t_72[k] = f_0 * pdd_18[k]
                  + f_4 * pc_z[k] * pfd_42[k];

        t_73[k] = f_0 * sfd_45[k]
                  + f_4 * pc_x[k] * pfd_45[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pa_x, pc_x, pc_y, pc_z, sff0_76, sfd_47, \
                         sff1_76, pdd_21, pdd_26, pfd_44, pfd_45, \
                         pfd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_6 * pdd_26[k]
                  + f_4 * pc_y[k] * pfd_44[k];

        t_75[k] = f_0 * sfd_47[k]
                  + f_4 * pc_x[k] * pfd_47[k];

        t_76[k] = pa_x[k] * sff0_76[k]
                  - f_5 * pc_x[k] * sff1_76[k];

        t_77[k] = f_0 * pdd_21[k]
                  + f_4 * pc_z[k] * pfd_45[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, pa_x, pb_y, pc_x, pc_y, sff0_79, sff1_79, \
                         pdf0_50, pdd_29, pdd_30, pdf1_50, pfd_47, \
                         pfd_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_6 * pdd_29[k]
                  + f_4 * pc_y[k] * pfd_47[k];

        t_79[k] = pa_x[k] * sff0_79[k]
                  - f_5 * pc_x[k] * sff1_79[k];

        t_80[k] = pb_y[k] * pdf0_50[k]
                  - f_5 * pc_y[k] * pdf1_50[k];

        t_81[k] = f_0 * pdd_30[k]
                  + f_4 * pc_y[k] * pfd_48[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pc_x, pc_y, pc_z, sfd_51, sfd_53, pdd_24, \
                         pdd_32, pfd_48, pfd_50, pfd_51, pfd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_6 * pdd_24[k]
                  + f_4 * pc_z[k] * pfd_48[k];

        t_83[k] = f_0 * sfd_51[k]
                  + f_4 * pc_x[k] * pfd_51[k];

        t_84[k] = f_0 * pdd_32[k]
                  + f_4 * pc_y[k] * pfd_50[k];

        t_85[k] = f_0 * sfd_53[k]
                  + f_4 * pc_x[k] * pfd_53[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pa_x, pc_x, pc_y, pc_z, sff0_86, sff0_89, \
                         sff1_86, sff1_89, pdd_27, pdd_35, pfd_51, \
                         pfd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = pa_x[k] * sff0_86[k]
                  - f_5 * pc_x[k] * sff1_86[k];

        t_87[k] = f_6 * pdd_27[k]
                  + f_4 * pc_z[k] * pfd_51[k];

        t_88[k] = f_0 * pdd_35[k]
                  + f_4 * pc_y[k] * pfd_53[k];

        t_89[k] = pa_x[k] * sff0_89[k]
                  - f_5 * pc_x[k] * sff1_89[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pa_x, pc_x, pc_y, pc_z, sff0_90, sfd_54, \
                         sfd_57, sff1_90, pdd_30, pfd_54, pfd_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = pa_x[k] * sff0_90[k]
                  + f_1 * sfd_54[k]
                  - f_5 * pc_x[k] * sff1_90[k];

        t_91[k] = f_4 * pc_y[k] * pfd_54[k];

        t_92[k] = f_1 * pdd_30[k]
                  + f_4 * pc_z[k] * pfd_54[k];

        t_93[k] = f_0 * sfd_57[k]
                  + f_4 * pc_x[k] * pfd_57[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, pa_x, pc_x, pc_y, pc_z, sff0_96, \
                         sfd_59, sff1_96, pdd_33, pfd_56, pfd_57, \
                         pfd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_4 * pc_y[k] * pfd_56[k];

        t_95[k] = f_0 * sfd_59[k]
                  + f_4 * pc_x[k] * pfd_59[k];

        t_96[k] = pa_x[k] * sff0_96[k]
                  - f_5 * pc_x[k] * sff1_96[k];

        t_97[k] = f_1 * pdd_33[k]
                  + f_4 * pc_z[k] * pfd_57[k];

        t_98[k] = f_4 * pc_y[k] * pfd_59[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pa_x, pa_y, pc_x, pc_y, sff0_0, sff0_1, sff0_99, \
                         sfd_0, sff1_0, sff1_1, sff1_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = pa_x[k] * sff0_99[k]
                  - f_5 * pc_x[k] * sff1_99[k];

        t_100[k] = pa_y[k] * sff0_0[k]
                   - f_5 * pc_y[k] * sff1_0[k];

        t_101[k] = pa_y[k] * sff0_1[k]
                   + f_0 * sfd_0[k]
                   - f_5 * pc_y[k] * sff1_1[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, pc_x, pc_z, pdd_39, pdd_40, pdd_41, \
                         pfd_60, pfd_63, pfd_64, pfd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_4 * pc_z[k] * pfd_60[k];

        t_103[k] = f_1 * pdd_39[k]
                   + f_4 * pc_x[k] * pfd_63[k];

        t_104[k] = f_1 * pdd_40[k]
                   + f_4 * pc_x[k] * pfd_64[k];

        t_105[k] = f_1 * pdd_41[k]
                   + f_4 * pc_x[k] * pfd_65[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pa_y, pc_y, pc_z, sff0_6, sff0_9, sfd_3, \
                         sfd_5, sff1_6, sff1_9, pfd_63, pfd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = pa_y[k] * sff0_6[k]
                   + f_1 * sfd_3[k]
                   - f_5 * pc_y[k] * sff1_6[k];

        t_107[k] = f_4 * pc_z[k] * pfd_63[k];

        t_108[k] = f_0 * sfd_5[k]
                   + f_4 * pc_y[k] * pfd_65[k];

        t_109[k] = pa_y[k] * sff0_9[k]
                   - f_5 * pc_y[k] * sff1_9[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pc_x, pc_z, pdd_42, pdd_43, pdd_45, \
                         pfp0_33, pfp0_34, pfp1_33, pfp1_34, pfd_66, pfd_67, \
                         pfd_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_6 * pdd_42[k]
                   + f_2 * pfp0_33[k]
                   - f_3 * pfp1_33[k]
                   + f_4 * pc_x[k] * pfd_66[k];

        t_111[k] = f_6 * pdd_43[k]
                   + f_7 * pfp0_34[k]
                   - f_8 * pfp1_34[k]
                   + f_4 * pc_x[k] * pfd_67[k];

        t_112[k] = f_4 * pc_z[k] * pfd_66[k];

        t_113[k] = f_6 * pdd_45[k]
                   + f_4 * pc_x[k] * pfd_69[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, t_117, pb_x, pc_x, pc_z, ppf0_46, ppf1_46, \
                         pdf0_76, pdd_46, pdd_47, pdf1_76, pfd_69, pfd_70, \
                         pfd_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_6 * pdd_46[k]
                   + f_4 * pc_x[k] * pfd_70[k];

        t_115[k] = f_6 * pdd_47[k]
                   + f_4 * pc_x[k] * pfd_71[k];

        t_116[k] = f_9 * ppf0_46[k]
                   - f_10 * ppf1_46[k]
                   + pb_x[k] * pdf0_76[k]
                   - f_5 * pc_x[k] * pdf1_76[k];

        t_117[k] = f_4 * pc_z[k] * pfd_69[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pa_y, pc_y, pc_z, sff0_20, sfd_11, sff1_20, \
                         pdd_41, pfp0_35, pfp1_35, pfd_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_0 * sfd_11[k]
                   + f_0 * pdd_41[k]
                   + f_4 * pc_y[k] * pfd_71[k];

        t_119[k] = f_2 * pfp0_35[k]
                   - f_3 * pfp1_35[k]
                   + f_4 * pc_z[k] * pfd_71[k];

        t_120[k] = pa_y[k] * sff0_20[k]
                   - f_5 * pc_y[k] * sff1_20[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pb_z, pc_x, pc_z, pdf0_61, pdd_36, \
                         pdd_51, pdd_52, pdf1_61, pfd_72, pfd_75, \
                         pfd_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = pb_z[k] * pdf0_61[k]
                   - f_5 * pc_z[k] * pdf1_61[k];

        t_122[k] = f_0 * pdd_36[k]
                   + f_4 * pc_z[k] * pfd_72[k];

        t_123[k] = f_6 * pdd_51[k]
                   + f_4 * pc_x[k] * pfd_75[k];

        t_124[k] = f_6 * pdd_52[k]
                   + f_4 * pc_x[k] * pfd_76[k];
    }

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pb_z, pc_x, pc_y, pc_z, sfd_17, pdf0_66, \
                         pdd_39, pdd_53, pdf1_66, pfd_75, pfd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_6 * pdd_53[k]
                   + f_4 * pc_x[k] * pfd_77[k];

        t_126[k] = pb_z[k] * pdf0_66[k]
                   - f_5 * pc_z[k] * pdf1_66[k];

        t_127[k] = f_0 * pdd_39[k]
                   + f_4 * pc_z[k] * pfd_75[k];

        t_128[k] = f_0 * sfd_17[k]
                   + f_4 * pc_y[k] * pfd_77[k];
    }
}

static auto
compute_prim_pff_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sff0, const size_t sfd,
                                                          const size_t sff1, const size_t ppf0,
                                                          const size_t ppf1, const size_t pdf0,
                                                          const size_t pdd, const size_t pdf1,
                                                          const size_t pfp0, const size_t pfp1,
                                                          const size_t pfd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.5 / q;
    const auto f_2 = 1.0 / gamma;
    const auto f_3 = p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = gamma / q;
    const auto f_6 = 1.0 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_9 = 0.5 / p;
    const auto f_10 = 0.5 * gamma / (p * q);

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
    auto *t_249 = buffer.data(target + 249);
    auto *t_250 = buffer.data(target + 250);
    auto *t_251 = buffer.data(target + 251);
    auto *t_252 = buffer.data(target + 252);
    auto *t_253 = buffer.data(target + 253);
    auto *t_254 = buffer.data(target + 254);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sff0_0 = buffer.data(sff0 + 0);
    const auto *sff0_2 = buffer.data(sff0 + 2);
    const auto *sff0_6 = buffer.data(sff0 + 6);
    const auto *sff0_9 = buffer.data(sff0 + 9);
    const auto *sff0_10 = buffer.data(sff0 + 10);
    const auto *sff0_16 = buffer.data(sff0 + 16);
    const auto *sff0_17 = buffer.data(sff0 + 17);
    const auto *sff0_29 = buffer.data(sff0 + 29);
    const auto *sff0_30 = buffer.data(sff0 + 30);
    const auto *sff0_32 = buffer.data(sff0 + 32);
    const auto *sff0_36 = buffer.data(sff0 + 36);
    const auto *sff0_50 = buffer.data(sff0 + 50);
    const auto *sff0_51 = buffer.data(sff0 + 51);
    const auto *sff0_59 = buffer.data(sff0 + 59);
    const auto *sff0_90 = buffer.data(sff0 + 90);
    const auto *sff0_96 = buffer.data(sff0 + 96);
    const auto *sff0_99 = buffer.data(sff0 + 99);

    const auto *sfd_0 = buffer.data(sfd + 0);
    const auto *sfd_5 = buffer.data(sfd + 5);
    const auto *sfd_9 = buffer.data(sfd + 9);
    const auto *sfd_18 = buffer.data(sfd + 18);
    const auto *sfd_30 = buffer.data(sfd + 30);
    const auto *sfd_35 = buffer.data(sfd + 35);
    const auto *sfd_39 = buffer.data(sfd + 39);
    const auto *sfd_41 = buffer.data(sfd + 41);
    const auto *sfd_47 = buffer.data(sfd + 47);
    const auto *sfd_51 = buffer.data(sfd + 51);
    const auto *sfd_53 = buffer.data(sfd + 53);
    const auto *sfd_57 = buffer.data(sfd + 57);
    const auto *sfd_59 = buffer.data(sfd + 59);

    const auto *sff1_0 = buffer.data(sff1 + 0);
    const auto *sff1_2 = buffer.data(sff1 + 2);
    const auto *sff1_6 = buffer.data(sff1 + 6);
    const auto *sff1_9 = buffer.data(sff1 + 9);
    const auto *sff1_10 = buffer.data(sff1 + 10);
    const auto *sff1_16 = buffer.data(sff1 + 16);
    const auto *sff1_17 = buffer.data(sff1 + 17);
    const auto *sff1_29 = buffer.data(sff1 + 29);
    const auto *sff1_30 = buffer.data(sff1 + 30);
    const auto *sff1_32 = buffer.data(sff1 + 32);
    const auto *sff1_36 = buffer.data(sff1 + 36);
    const auto *sff1_50 = buffer.data(sff1 + 50);
    const auto *sff1_51 = buffer.data(sff1 + 51);
    const auto *sff1_59 = buffer.data(sff1 + 59);
    const auto *sff1_90 = buffer.data(sff1 + 90);
    const auto *sff1_96 = buffer.data(sff1 + 96);
    const auto *sff1_99 = buffer.data(sff1 + 99);

    const auto *ppf0_89 = buffer.data(ppf0 + 89);

    const auto *ppf1_89 = buffer.data(ppf1 + 89);

    const auto *pdf0_71 = buffer.data(pdf0 + 71);
    const auto *pdf0_90 = buffer.data(pdf0 + 90);
    const auto *pdf0_91 = buffer.data(pdf0 + 91);
    const auto *pdf0_96 = buffer.data(pdf0 + 96);
    const auto *pdf0_98 = buffer.data(pdf0 + 98);
    const auto *pdf0_99 = buffer.data(pdf0 + 99);
    const auto *pdf0_106 = buffer.data(pdf0 + 106);
    const auto *pdf0_108 = buffer.data(pdf0 + 108);
    const auto *pdf0_109 = buffer.data(pdf0 + 109);
    const auto *pdf0_116 = buffer.data(pdf0 + 116);
    const auto *pdf0_122 = buffer.data(pdf0 + 122);
    const auto *pdf0_129 = buffer.data(pdf0 + 129);
    const auto *pdf0_140 = buffer.data(pdf0 + 140);
    const auto *pdf0_142 = buffer.data(pdf0 + 142);
    const auto *pdf0_149 = buffer.data(pdf0 + 149);
    const auto *pdf0_157 = buffer.data(pdf0 + 157);
    const auto *pdf0_159 = buffer.data(pdf0 + 159);
    const auto *pdf0_166 = buffer.data(pdf0 + 166);
    const auto *pdf0_167 = buffer.data(pdf0 + 167);
    const auto *pdf0_169 = buffer.data(pdf0 + 169);
    const auto *pdf0_170 = buffer.data(pdf0 + 170);
    const auto *pdf0_172 = buffer.data(pdf0 + 172);

    const auto *pdd_42 = buffer.data(pdd + 42);
    const auto *pdd_45 = buffer.data(pdd + 45);
    const auto *pdd_48 = buffer.data(pdd + 48);
    const auto *pdd_51 = buffer.data(pdd + 51);
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
    const auto *pdd_69 = buffer.data(pdd + 69);
    const auto *pdd_70 = buffer.data(pdd + 70);
    const auto *pdd_71 = buffer.data(pdd + 71);
    const auto *pdd_72 = buffer.data(pdd + 72);
    const auto *pdd_75 = buffer.data(pdd + 75);
    const auto *pdd_76 = buffer.data(pdd + 76);
    const auto *pdd_77 = buffer.data(pdd + 77);
    const auto *pdd_78 = buffer.data(pdd + 78);
    const auto *pdd_81 = buffer.data(pdd + 81);
    const auto *pdd_82 = buffer.data(pdd + 82);
    const auto *pdd_83 = buffer.data(pdd + 83);
    const auto *pdd_84 = buffer.data(pdd + 84);
    const auto *pdd_86 = buffer.data(pdd + 86);
    const auto *pdd_87 = buffer.data(pdd + 87);
    const auto *pdd_88 = buffer.data(pdd + 88);
    const auto *pdd_89 = buffer.data(pdd + 89);
    const auto *pdd_93 = buffer.data(pdd + 93);
    const auto *pdd_94 = buffer.data(pdd + 94);
    const auto *pdd_95 = buffer.data(pdd + 95);
    const auto *pdd_99 = buffer.data(pdd + 99);
    const auto *pdd_100 = buffer.data(pdd + 100);
    const auto *pdd_101 = buffer.data(pdd + 101);
    const auto *pdd_102 = buffer.data(pdd + 102);
    const auto *pdd_104 = buffer.data(pdd + 104);
    const auto *pdd_105 = buffer.data(pdd + 105);
    const auto *pdd_106 = buffer.data(pdd + 106);

    const auto *pdf1_71 = buffer.data(pdf1 + 71);
    const auto *pdf1_90 = buffer.data(pdf1 + 90);
    const auto *pdf1_91 = buffer.data(pdf1 + 91);
    const auto *pdf1_96 = buffer.data(pdf1 + 96);
    const auto *pdf1_98 = buffer.data(pdf1 + 98);
    const auto *pdf1_99 = buffer.data(pdf1 + 99);
    const auto *pdf1_106 = buffer.data(pdf1 + 106);
    const auto *pdf1_108 = buffer.data(pdf1 + 108);
    const auto *pdf1_109 = buffer.data(pdf1 + 109);
    const auto *pdf1_116 = buffer.data(pdf1 + 116);
    const auto *pdf1_122 = buffer.data(pdf1 + 122);
    const auto *pdf1_129 = buffer.data(pdf1 + 129);
    const auto *pdf1_140 = buffer.data(pdf1 + 140);
    const auto *pdf1_142 = buffer.data(pdf1 + 142);
    const auto *pdf1_149 = buffer.data(pdf1 + 149);
    const auto *pdf1_157 = buffer.data(pdf1 + 157);
    const auto *pdf1_159 = buffer.data(pdf1 + 159);
    const auto *pdf1_166 = buffer.data(pdf1 + 166);
    const auto *pdf1_167 = buffer.data(pdf1 + 167);
    const auto *pdf1_169 = buffer.data(pdf1 + 169);
    const auto *pdf1_170 = buffer.data(pdf1 + 170);
    const auto *pdf1_172 = buffer.data(pdf1 + 172);

    const auto *pfp0_42 = buffer.data(pfp0 + 42);
    const auto *pfp0_48 = buffer.data(pfp0 + 48);
    const auto *pfp0_49 = buffer.data(pfp0 + 49);
    const auto *pfp0_50 = buffer.data(pfp0 + 50);
    const auto *pfp0_53 = buffer.data(pfp0 + 53);
    const auto *pfp0_54 = buffer.data(pfp0 + 54);
    const auto *pfp0_55 = buffer.data(pfp0 + 55);
    const auto *pfp0_56 = buffer.data(pfp0 + 56);
    const auto *pfp0_58 = buffer.data(pfp0 + 58);
    const auto *pfp0_62 = buffer.data(pfp0 + 62);
    const auto *pfp0_66 = buffer.data(pfp0 + 66);
    const auto *pfp0_67 = buffer.data(pfp0 + 67);
    const auto *pfp0_68 = buffer.data(pfp0 + 68);

    const auto *pfp1_42 = buffer.data(pfp1 + 42);
    const auto *pfp1_48 = buffer.data(pfp1 + 48);
    const auto *pfp1_49 = buffer.data(pfp1 + 49);
    const auto *pfp1_50 = buffer.data(pfp1 + 50);
    const auto *pfp1_53 = buffer.data(pfp1 + 53);
    const auto *pfp1_54 = buffer.data(pfp1 + 54);
    const auto *pfp1_55 = buffer.data(pfp1 + 55);
    const auto *pfp1_56 = buffer.data(pfp1 + 56);
    const auto *pfp1_58 = buffer.data(pfp1 + 58);
    const auto *pfp1_62 = buffer.data(pfp1 + 62);
    const auto *pfp1_66 = buffer.data(pfp1 + 66);
    const auto *pfp1_67 = buffer.data(pfp1 + 67);
    const auto *pfp1_68 = buffer.data(pfp1 + 68);

    const auto *pfd_78 = buffer.data(pfd + 78);
    const auto *pfd_81 = buffer.data(pfd + 81);
    const auto *pfd_82 = buffer.data(pfd + 82);
    const auto *pfd_83 = buffer.data(pfd + 83);
    const auto *pfd_84 = buffer.data(pfd + 84);
    const auto *pfd_87 = buffer.data(pfd + 87);
    const auto *pfd_88 = buffer.data(pfd + 88);
    const auto *pfd_89 = buffer.data(pfd + 89);
    const auto *pfd_90 = buffer.data(pfd + 90);
    const auto *pfd_93 = buffer.data(pfd + 93);
    const auto *pfd_94 = buffer.data(pfd + 94);
    const auto *pfd_95 = buffer.data(pfd + 95);
    const auto *pfd_96 = buffer.data(pfd + 96);
    const auto *pfd_97 = buffer.data(pfd + 97);
    const auto *pfd_99 = buffer.data(pfd + 99);
    const auto *pfd_100 = buffer.data(pfd + 100);
    const auto *pfd_101 = buffer.data(pfd + 101);
    const auto *pfd_102 = buffer.data(pfd + 102);
    const auto *pfd_105 = buffer.data(pfd + 105);
    const auto *pfd_106 = buffer.data(pfd + 106);
    const auto *pfd_107 = buffer.data(pfd + 107);
    const auto *pfd_108 = buffer.data(pfd + 108);
    const auto *pfd_109 = buffer.data(pfd + 109);
    const auto *pfd_111 = buffer.data(pfd + 111);
    const auto *pfd_112 = buffer.data(pfd + 112);
    const auto *pfd_113 = buffer.data(pfd + 113);
    const auto *pfd_114 = buffer.data(pfd + 114);
    const auto *pfd_115 = buffer.data(pfd + 115);
    const auto *pfd_117 = buffer.data(pfd + 117);
    const auto *pfd_118 = buffer.data(pfd + 118);
    const auto *pfd_119 = buffer.data(pfd + 119);
    const auto *pfd_120 = buffer.data(pfd + 120);
    const auto *pfd_123 = buffer.data(pfd + 123);
    const auto *pfd_124 = buffer.data(pfd + 124);
    const auto *pfd_125 = buffer.data(pfd + 125);
    const auto *pfd_126 = buffer.data(pfd + 126);
    const auto *pfd_129 = buffer.data(pfd + 129);
    const auto *pfd_130 = buffer.data(pfd + 130);
    const auto *pfd_131 = buffer.data(pfd + 131);
    const auto *pfd_132 = buffer.data(pfd + 132);
    const auto *pfd_134 = buffer.data(pfd + 134);
    const auto *pfd_135 = buffer.data(pfd + 135);
    const auto *pfd_136 = buffer.data(pfd + 136);
    const auto *pfd_137 = buffer.data(pfd + 137);
    const auto *pfd_138 = buffer.data(pfd + 138);
    const auto *pfd_141 = buffer.data(pfd + 141);
    const auto *pfd_142 = buffer.data(pfd + 142);
    const auto *pfd_143 = buffer.data(pfd + 143);
    const auto *pfd_144 = buffer.data(pfd + 144);
    const auto *pfd_147 = buffer.data(pfd + 147);
    const auto *pfd_148 = buffer.data(pfd + 148);
    const auto *pfd_149 = buffer.data(pfd + 149);
    const auto *pfd_150 = buffer.data(pfd + 150);
    const auto *pfd_153 = buffer.data(pfd + 153);
    const auto *pfd_154 = buffer.data(pfd + 154);

#pragma omp simd aligned(t_129, t_130, t_131, pa_y, pb_x, pc_x, pc_y, sff0_29, sff1_29, \
                         pdf0_90, pdf0_91, pdd_54, pdd_55, pdf1_90, \
                         pdf1_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = pa_y[k] * sff0_29[k]
                   - f_5 * pc_y[k] * sff1_29[k];

        t_130[k] = pb_x[k] * pdf0_90[k]
                   + f_1 * pdd_54[k]
                   - f_5 * pc_x[k] * pdf1_90[k];

        t_131[k] = pb_x[k] * pdf0_91[k]
                   + f_6 * pdd_55[k]
                   - f_5 * pc_x[k] * pdf1_91[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pc_x, pc_z, pdd_57, pdd_58, pdd_59, \
                         pfd_78, pfd_81, pfd_82, pfd_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_4 * pc_z[k] * pfd_78[k];

        t_133[k] = f_0 * pdd_57[k]
                   + f_4 * pc_x[k] * pfd_81[k];

        t_134[k] = f_0 * pdd_58[k]
                   + f_4 * pc_x[k] * pfd_82[k];

        t_135[k] = f_0 * pdd_59[k]
                   + f_4 * pc_x[k] * pfd_83[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pb_x, pc_x, pc_z, pdf0_96, pdf0_98, \
                         pdf0_99, pdf1_96, pdf1_98, pdf1_99, pfd_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = pb_x[k] * pdf0_96[k]
                   - f_5 * pc_x[k] * pdf1_96[k];

        t_137[k] = f_4 * pc_z[k] * pfd_81[k];

        t_138[k] = pb_x[k] * pdf0_98[k]
                   - f_5 * pc_x[k] * pdf1_98[k];

        t_139[k] = pb_x[k] * pdf0_99[k]
                   - f_5 * pc_x[k] * pdf1_99[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pb_z, pc_x, pc_z, pdf0_71, pdd_42, \
                         pdd_60, pdd_63, pdf1_71, pfp0_42, pfp1_42, pfd_84, \
                         pfd_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_0 * pdd_60[k]
                   + f_2 * pfp0_42[k]
                   - f_3 * pfp1_42[k]
                   + f_4 * pc_x[k] * pfd_84[k];

        t_141[k] = pb_z[k] * pdf0_71[k]
                   - f_5 * pc_z[k] * pdf1_71[k];

        t_142[k] = f_0 * pdd_42[k]
                   + f_4 * pc_z[k] * pfd_84[k];

        t_143[k] = f_0 * pdd_63[k]
                   + f_4 * pc_x[k] * pfd_87[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pb_x, pc_x, pc_z, pdf0_106, pdd_45, \
                         pdd_64, pdd_65, pdf1_106, pfd_87, pfd_88, \
                         pfd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_0 * pdd_64[k]
                   + f_4 * pc_x[k] * pfd_88[k];

        t_145[k] = f_0 * pdd_65[k]
                   + f_4 * pc_x[k] * pfd_89[k];

        t_146[k] = pb_x[k] * pdf0_106[k]
                   - f_5 * pc_x[k] * pdf1_106[k];

        t_147[k] = f_0 * pdd_45[k]
                   + f_4 * pc_z[k] * pfd_87[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, pa_y, pb_x, pc_x, pc_y, sff0_50, sff1_50, \
                         pdf0_108, pdf0_109, pdf1_108, pdf1_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = pb_x[k] * pdf0_108[k]
                   - f_5 * pc_x[k] * pdf1_108[k];

        t_149[k] = pb_x[k] * pdf0_109[k]
                   - f_5 * pc_x[k] * pdf1_109[k];

        t_150[k] = pa_y[k] * sff0_50[k]
                   - f_5 * pc_y[k] * sff1_50[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, pa_y, pc_x, pc_y, pc_z, sff0_51, sfd_30, \
                         sff1_51, pdd_48, pdd_69, pfd_90, pfd_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = pa_y[k] * sff0_51[k]
                   + f_0 * sfd_30[k]
                   - f_5 * pc_y[k] * sff1_51[k];

        t_152[k] = f_6 * pdd_48[k]
                   + f_4 * pc_z[k] * pfd_90[k];

        t_153[k] = f_0 * pdd_69[k]
                   + f_4 * pc_x[k] * pfd_93[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pb_x, pc_x, pc_z, pdf0_116, pdd_51, \
                         pdd_70, pdd_71, pdf1_116, pfd_93, pfd_94, \
                         pfd_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_0 * pdd_70[k]
                   + f_4 * pc_x[k] * pfd_94[k];

        t_155[k] = f_0 * pdd_71[k]
                   + f_4 * pc_x[k] * pfd_95[k];

        t_156[k] = pb_x[k] * pdf0_116[k]
                   - f_5 * pc_x[k] * pdf1_116[k];

        t_157[k] = f_6 * pdd_51[k]
                   + f_4 * pc_z[k] * pfd_93[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, pa_y, pc_x, pc_y, sff0_59, sfd_35, sff1_59, \
                         pfp0_48, pfp1_48, pfd_95, pfd_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_0 * sfd_35[k]
                   + f_4 * pc_y[k] * pfd_95[k];

        t_159[k] = pa_y[k] * sff0_59[k]
                   - f_5 * pc_y[k] * sff1_59[k];

        t_160[k] = f_2 * pfp0_48[k]
                   - f_3 * pfp1_48[k]
                   + f_4 * pc_x[k] * pfd_96[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, t_165, pc_x, pc_z, pfp0_49, pfp1_49, \
                         pfd_96, pfd_97, pfd_99, pfd_100, pfd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_7 * pfp0_49[k]
                   - f_8 * pfp1_49[k]
                   + f_4 * pc_x[k] * pfd_97[k];

        t_162[k] = f_4 * pc_z[k] * pfd_96[k];

        t_163[k] = f_4 * pc_x[k] * pfd_99[k];

        t_164[k] = f_4 * pc_x[k] * pfd_100[k];

        t_165[k] = f_4 * pc_x[k] * pfd_101[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pc_y, pc_z, sfd_39, sfd_41, pdd_57, \
                         pdd_59, pfp0_49, pfp0_50, pfp1_49, pfp1_50, pfd_99, \
                         pfd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_0 * sfd_39[k]
                   + f_1 * pdd_57[k]
                   + f_2 * pfp0_49[k]
                   - f_3 * pfp1_49[k]
                   + f_4 * pc_y[k] * pfd_99[k];

        t_167[k] = f_4 * pc_z[k] * pfd_99[k];

        t_168[k] = f_0 * sfd_41[k]
                   + f_1 * pdd_59[k]
                   + f_4 * pc_y[k] * pfd_101[k];

        t_169[k] = f_2 * pfp0_50[k]
                   - f_3 * pfp1_50[k]
                   + f_4 * pc_z[k] * pfd_101[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, t_174, pb_z, pc_x, pc_z, pdf0_90, \
                         pdf0_91, pdd_54, pdf1_90, pdf1_91, pfd_102, pfd_105, \
                         pfd_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = pb_z[k] * pdf0_90[k]
                   - f_5 * pc_z[k] * pdf1_90[k];

        t_171[k] = pb_z[k] * pdf0_91[k]
                   - f_5 * pc_z[k] * pdf1_91[k];

        t_172[k] = f_0 * pdd_54[k]
                   + f_4 * pc_z[k] * pfd_102[k];

        t_173[k] = f_4 * pc_x[k] * pfd_105[k];

        t_174[k] = f_4 * pc_x[k] * pfd_106[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, t_178, pb_z, pc_x, pc_y, pc_z, sfd_47, pdf0_96, \
                         pdd_57, pdd_65, pdf1_96, pfd_105, pfd_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_4 * pc_x[k] * pfd_107[k];

        t_176[k] = pb_z[k] * pdf0_96[k]
                   - f_5 * pc_z[k] * pdf1_96[k];

        t_177[k] = f_0 * pdd_57[k]
                   + f_4 * pc_z[k] * pfd_105[k];

        t_178[k] = f_0 * sfd_47[k]
                   + f_6 * pdd_65[k]
                   + f_4 * pc_y[k] * pfd_107[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, pc_x, pc_z, pdd_59, pfp0_53, pfp0_54, pfp0_55, \
                         pfp1_53, pfp1_54, pfp1_55, pfd_107, pfd_108, \
                         pfd_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_0 * pdd_59[k]
                   + f_2 * pfp0_53[k]
                   - f_3 * pfp1_53[k]
                   + f_4 * pc_z[k] * pfd_107[k];

        t_180[k] = f_2 * pfp0_54[k]
                   - f_3 * pfp1_54[k]
                   + f_4 * pc_x[k] * pfd_108[k];

        t_181[k] = f_7 * pfp0_55[k]
                   - f_8 * pfp1_55[k]
                   + f_4 * pc_x[k] * pfd_109[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, t_185, pc_x, pc_z, pdd_60, pfd_108, pfd_111, \
                         pfd_112, pfd_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = f_6 * pdd_60[k]
                   + f_4 * pc_z[k] * pfd_108[k];

        t_183[k] = f_4 * pc_x[k] * pfd_111[k];

        t_184[k] = f_4 * pc_x[k] * pfd_112[k];

        t_185[k] = f_4 * pc_x[k] * pfd_113[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, pc_y, pc_z, sfd_51, sfd_53, pdd_63, pdd_69, \
                         pdd_71, pfp0_55, pfp1_55, pfd_111, pfd_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_0 * sfd_51[k]
                   + f_0 * pdd_69[k]
                   + f_2 * pfp0_55[k]
                   - f_3 * pfp1_55[k]
                   + f_4 * pc_y[k] * pfd_111[k];

        t_187[k] = f_6 * pdd_63[k]
                   + f_4 * pc_z[k] * pfd_111[k];

        t_188[k] = f_0 * sfd_53[k]
                   + f_0 * pdd_71[k]
                   + f_4 * pc_y[k] * pfd_113[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pa_y, pc_x, pc_y, pc_z, sff0_90, sff1_90, \
                         pdd_65, pfp0_56, pfp0_58, pfp1_56, pfp1_58, pfd_113, \
                         pfd_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_6 * pdd_65[k]
                   + f_2 * pfp0_56[k]
                   - f_3 * pfp1_56[k]
                   + f_4 * pc_z[k] * pfd_113[k];

        t_190[k] = pa_y[k] * sff0_90[k]
                   - f_5 * pc_y[k] * sff1_90[k];

        t_191[k] = f_7 * pfp0_58[k]
                   - f_8 * pfp1_58[k]
                   + f_4 * pc_x[k] * pfd_115[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pc_x, pc_z, pdd_66, pfd_114, pfd_117, \
                         pfd_118, pfd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_1 * pdd_66[k]
                   + f_4 * pc_z[k] * pfd_114[k];

        t_193[k] = f_4 * pc_x[k] * pfd_117[k];

        t_194[k] = f_4 * pc_x[k] * pfd_118[k];

        t_195[k] = f_4 * pc_x[k] * pfd_119[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pa_y, pc_y, pc_z, sff0_96, sff0_99, \
                         sfd_57, sfd_59, sff1_96, sff1_99, pdd_69, pfd_117, \
                         pfd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = pa_y[k] * sff0_96[k]
                   + f_1 * sfd_57[k]
                   - f_5 * pc_y[k] * sff1_96[k];

        t_197[k] = f_1 * pdd_69[k]
                   + f_4 * pc_z[k] * pfd_117[k];

        t_198[k] = f_0 * sfd_59[k]
                   + f_4 * pc_y[k] * pfd_119[k];

        t_199[k] = pa_y[k] * sff0_99[k]
                   - f_5 * pc_y[k] * sff1_99[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, t_203, pa_z, pc_x, pc_y, pc_z, sff0_0, sff0_2, \
                         sfd_0, sff1_0, sff1_2, pdd_75, pfd_120, \
                         pfd_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = pa_z[k] * sff0_0[k]
                   - f_5 * pc_z[k] * sff1_0[k];

        t_201[k] = f_4 * pc_y[k] * pfd_120[k];

        t_202[k] = pa_z[k] * sff0_2[k]
                   + f_0 * sfd_0[k]
                   - f_5 * pc_z[k] * sff1_2[k];

        t_203[k] = f_1 * pdd_75[k]
                   + f_4 * pc_x[k] * pfd_123[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pa_z, pc_x, pc_y, pc_z, sff0_6, sff1_6, \
                         pdd_76, pdd_77, pfp0_62, pfp1_62, pfd_124, \
                         pfd_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = f_1 * pdd_76[k]
                   + f_4 * pc_x[k] * pfd_124[k];

        t_205[k] = f_1 * pdd_77[k]
                   + f_4 * pc_x[k] * pfd_125[k];

        t_206[k] = pa_z[k] * sff0_6[k]
                   - f_5 * pc_z[k] * sff1_6[k];

        t_207[k] = f_7 * pfp0_62[k]
                   - f_8 * pfp1_62[k]
                   + f_4 * pc_y[k] * pfd_124[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pa_z, pc_y, pc_z, sff0_9, sff0_10, sfd_5, \
                         sff1_9, sff1_10, pdd_72, pfd_125, pfd_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_4 * pc_y[k] * pfd_125[k];

        t_209[k] = pa_z[k] * sff0_9[k]
                   + f_1 * sfd_5[k]
                   - f_5 * pc_z[k] * sff1_9[k];

        t_210[k] = pa_z[k] * sff0_10[k]
                   - f_5 * pc_z[k] * sff1_10[k];

        t_211[k] = f_0 * pdd_72[k]
                   + f_4 * pc_y[k] * pfd_126[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, t_215, pb_y, pc_x, pc_y, pdf0_122, pdd_81, \
                         pdd_82, pdd_83, pdf1_122, pfd_129, pfd_130, \
                         pfd_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = pb_y[k] * pdf0_122[k]
                   - f_5 * pc_y[k] * pdf1_122[k];

        t_213[k] = f_6 * pdd_81[k]
                   + f_4 * pc_x[k] * pfd_129[k];

        t_214[k] = f_6 * pdd_82[k]
                   + f_4 * pc_x[k] * pfd_130[k];

        t_215[k] = f_6 * pdd_83[k]
                   + f_4 * pc_x[k] * pfd_131[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, pa_z, pc_y, pc_z, sff0_16, sff0_17, sfd_9, \
                         sff1_16, sff1_17, pdd_77, pfd_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = pa_z[k] * sff0_16[k]
                   - f_5 * pc_z[k] * sff1_16[k];

        t_217[k] = pa_z[k] * sff0_17[k]
                   + f_0 * sfd_9[k]
                   - f_5 * pc_z[k] * sff1_17[k];

        t_218[k] = f_0 * pdd_77[k]
                   + f_4 * pc_y[k] * pfd_131[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, pb_y, pc_x, pc_y, pdf0_129, pdd_84, pdf1_129, \
                         pfp0_66, pfp1_66, pfd_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = pb_y[k] * pdf0_129[k]
                   - f_5 * pc_y[k] * pdf1_129[k];

        t_220[k] = f_6 * pdd_84[k]
                   + f_2 * pfp0_66[k]
                   - f_3 * pfp1_66[k]
                   + f_4 * pc_x[k] * pfd_132[k];

        t_221[k] = f_4 * pc_y[k] * pfd_132[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pc_x, pdd_86, pdd_87, pdd_88, pdd_89, \
                         pfp0_68, pfp1_68, pfd_134, pfd_135, pfd_136, \
                         pfd_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = f_6 * pdd_86[k]
                   + f_7 * pfp0_68[k]
                   - f_8 * pfp1_68[k]
                   + f_4 * pc_x[k] * pfd_134[k];

        t_223[k] = f_6 * pdd_87[k]
                   + f_4 * pc_x[k] * pfd_135[k];

        t_224[k] = f_6 * pdd_88[k]
                   + f_4 * pc_x[k] * pfd_136[k];

        t_225[k] = f_6 * pdd_89[k]
                   + f_4 * pc_x[k] * pfd_137[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, pc_y, pfp0_67, pfp0_68, pfp1_67, pfp1_68, \
                         pfd_135, pfd_136, pfd_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = f_2 * pfp0_67[k]
                   - f_3 * pfp1_67[k]
                   + f_4 * pc_y[k] * pfd_135[k];

        t_227[k] = f_7 * pfp0_68[k]
                   - f_8 * pfp1_68[k]
                   + f_4 * pc_y[k] * pfd_136[k];

        t_228[k] = f_4 * pc_y[k] * pfd_137[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, pa_z, pb_x, pc_x, pc_y, pc_z, sff0_30, sff1_30, \
                         ppf0_89, ppf1_89, pdf0_149, pdd_78, pdf1_149, \
                         pfd_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_9 * ppf0_89[k]
                   - f_10 * ppf1_89[k]
                   + pb_x[k] * pdf0_149[k]
                   - f_5 * pc_x[k] * pdf1_149[k];

        t_230[k] = pa_z[k] * sff0_30[k]
                   - f_5 * pc_z[k] * sff1_30[k];

        t_231[k] = f_6 * pdd_78[k]
                   + f_4 * pc_y[k] * pfd_138[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, t_235, pa_z, pc_x, pc_z, sff0_32, sfd_18, \
                         sff1_32, pdd_93, pdd_94, pdd_95, pfd_141, pfd_142, \
                         pfd_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = pa_z[k] * sff0_32[k]
                   + f_0 * sfd_18[k]
                   - f_5 * pc_z[k] * sff1_32[k];

        t_233[k] = f_0 * pdd_93[k]
                   + f_4 * pc_x[k] * pfd_141[k];

        t_234[k] = f_0 * pdd_94[k]
                   + f_4 * pc_x[k] * pfd_142[k];

        t_235[k] = f_0 * pdd_95[k]
                   + f_4 * pc_x[k] * pfd_143[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, pa_z, pb_x, pc_x, pc_y, pc_z, sff0_36, sff1_36, \
                         pdf0_157, pdd_83, pdf1_157, pfd_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = pa_z[k] * sff0_36[k]
                   - f_5 * pc_z[k] * sff1_36[k];

        t_237[k] = pb_x[k] * pdf0_157[k]
                   - f_5 * pc_x[k] * pdf1_157[k];

        t_238[k] = f_6 * pdd_83[k]
                   + f_4 * pc_y[k] * pfd_143[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, pb_x, pb_y, pc_x, pc_y, pdf0_140, \
                         pdf0_142, pdf0_159, pdd_84, pdf1_140, pdf1_142, pdf1_159, \
                         pfd_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = pb_x[k] * pdf0_159[k]
                   - f_5 * pc_x[k] * pdf1_159[k];

        t_240[k] = pb_y[k] * pdf0_140[k]
                   - f_5 * pc_y[k] * pdf1_140[k];

        t_241[k] = f_0 * pdd_84[k]
                   + f_4 * pc_y[k] * pfd_144[k];

        t_242[k] = pb_y[k] * pdf0_142[k]
                   - f_5 * pc_y[k] * pdf1_142[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, t_246, pb_x, pc_x, pdf0_166, pdd_99, pdd_100, \
                         pdd_101, pdf1_166, pfd_147, pfd_148, pfd_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_0 * pdd_99[k]
                   + f_4 * pc_x[k] * pfd_147[k];

        t_244[k] = f_0 * pdd_100[k]
                   + f_4 * pc_x[k] * pfd_148[k];

        t_245[k] = f_0 * pdd_101[k]
                   + f_4 * pc_x[k] * pfd_149[k];

        t_246[k] = pb_x[k] * pdf0_166[k]
                   - f_5 * pc_x[k] * pdf1_166[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, pb_x, pc_x, pc_y, pdf0_167, pdf0_169, \
                         pdf0_170, pdd_89, pdd_102, pdf1_167, pdf1_169, pdf1_170, \
                         pfd_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = pb_x[k] * pdf0_167[k]
                   - f_5 * pc_x[k] * pdf1_167[k];

        t_248[k] = f_0 * pdd_89[k]
                   + f_4 * pc_y[k] * pfd_149[k];

        t_249[k] = pb_x[k] * pdf0_169[k]
                   - f_5 * pc_x[k] * pdf1_169[k];

        t_250[k] = pb_x[k] * pdf0_170[k]
                   + f_1 * pdd_102[k]
                   - f_5 * pc_x[k] * pdf1_170[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, t_254, pb_x, pc_x, pc_y, pdf0_172, pdd_104, \
                         pdd_105, pdd_106, pdf1_172, pfd_150, pfd_153, \
                         pfd_154 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_4 * pc_y[k] * pfd_150[k];

        t_252[k] = pb_x[k] * pdf0_172[k]
                   + f_6 * pdd_104[k]
                   - f_5 * pc_x[k] * pdf1_172[k];

        t_253[k] = f_0 * pdd_105[k]
                   + f_4 * pc_x[k] * pfd_153[k];

        t_254[k] = f_0 * pdd_106[k]
                   + f_4 * pc_x[k] * pfd_154[k];
    }
}

static auto
compute_prim_pff_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sff0, const size_t sfd,
                                                          const size_t sff1, const size_t ppf0,
                                                          const size_t ppf1, const size_t pdf0,
                                                          const size_t pdd, const size_t pdf1,
                                                          const size_t pfp0, const size_t pfp1,
                                                          const size_t pfd, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.5 / q;
    const auto f_2 = 1.0 / gamma;
    const auto f_3 = p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = gamma / q;
    const auto f_6 = 1.0 / q;
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);
    const auto f_9 = 0.5 / p;
    const auto f_10 = 0.5 * gamma / (p * q);

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

    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sff0_60 = buffer.data(sff0 + 60);
    const auto *sff0_66 = buffer.data(sff0 + 66);
    const auto *sff0_67 = buffer.data(sff0 + 67);
    const auto *sff0_69 = buffer.data(sff0 + 69);

    const auto *sfd_39 = buffer.data(sfd + 39);
    const auto *sfd_41 = buffer.data(sfd + 41);
    const auto *sfd_59 = buffer.data(sfd + 59);

    const auto *sff1_60 = buffer.data(sff1 + 60);
    const auto *sff1_66 = buffer.data(sff1 + 66);
    const auto *sff1_67 = buffer.data(sff1 + 67);
    const auto *sff1_69 = buffer.data(sff1 + 69);

    const auto *ppf0_89 = buffer.data(ppf0 + 89);

    const auto *ppf1_89 = buffer.data(ppf1 + 89);

    const auto *pdf0_169 = buffer.data(pdf0 + 169);
    const auto *pdf0_170 = buffer.data(pdf0 + 170);
    const auto *pdf0_172 = buffer.data(pdf0 + 172);
    const auto *pdf0_176 = buffer.data(pdf0 + 176);
    const auto *pdf0_177 = buffer.data(pdf0 + 177);
    const auto *pdf0_179 = buffer.data(pdf0 + 179);

    const auto *pdd_90 = buffer.data(pdd + 90);
    const auto *pdd_95 = buffer.data(pdd + 95);
    const auto *pdd_96 = buffer.data(pdd + 96);
    const auto *pdd_99 = buffer.data(pdd + 99);
    const auto *pdd_100 = buffer.data(pdd + 100);
    const auto *pdd_101 = buffer.data(pdd + 101);
    const auto *pdd_102 = buffer.data(pdd + 102);
    const auto *pdd_105 = buffer.data(pdd + 105);
    const auto *pdd_106 = buffer.data(pdd + 106);
    const auto *pdd_107 = buffer.data(pdd + 107);

    const auto *pdf1_169 = buffer.data(pdf1 + 169);
    const auto *pdf1_170 = buffer.data(pdf1 + 170);
    const auto *pdf1_172 = buffer.data(pdf1 + 172);
    const auto *pdf1_176 = buffer.data(pdf1 + 176);
    const auto *pdf1_177 = buffer.data(pdf1 + 177);
    const auto *pdf1_179 = buffer.data(pdf1 + 179);

    const auto *pfp0_80 = buffer.data(pfp0 + 80);
    const auto *pfp0_81 = buffer.data(pfp0 + 81);
    const auto *pfp0_82 = buffer.data(pfp0 + 82);
    const auto *pfp0_83 = buffer.data(pfp0 + 83);
    const auto *pfp0_87 = buffer.data(pfp0 + 87);
    const auto *pfp0_88 = buffer.data(pfp0 + 88);
    const auto *pfp0_89 = buffer.data(pfp0 + 89);

    const auto *pfp1_80 = buffer.data(pfp1 + 80);
    const auto *pfp1_81 = buffer.data(pfp1 + 81);
    const auto *pfp1_82 = buffer.data(pfp1 + 82);
    const auto *pfp1_83 = buffer.data(pfp1 + 83);
    const auto *pfp1_87 = buffer.data(pfp1 + 87);
    const auto *pfp1_88 = buffer.data(pfp1 + 88);
    const auto *pfp1_89 = buffer.data(pfp1 + 89);

    const auto *pfd_155 = buffer.data(pfd + 155);
    const auto *pfd_156 = buffer.data(pfd + 156);
    const auto *pfd_158 = buffer.data(pfd + 158);
    const auto *pfd_159 = buffer.data(pfd + 159);
    const auto *pfd_160 = buffer.data(pfd + 160);
    const auto *pfd_161 = buffer.data(pfd + 161);
    const auto *pfd_162 = buffer.data(pfd + 162);
    const auto *pfd_164 = buffer.data(pfd + 164);
    const auto *pfd_165 = buffer.data(pfd + 165);
    const auto *pfd_166 = buffer.data(pfd + 166);
    const auto *pfd_167 = buffer.data(pfd + 167);
    const auto *pfd_168 = buffer.data(pfd + 168);
    const auto *pfd_171 = buffer.data(pfd + 171);
    const auto *pfd_172 = buffer.data(pfd + 172);
    const auto *pfd_173 = buffer.data(pfd + 173);
    const auto *pfd_174 = buffer.data(pfd + 174);
    const auto *pfd_176 = buffer.data(pfd + 176);
    const auto *pfd_177 = buffer.data(pfd + 177);
    const auto *pfd_178 = buffer.data(pfd + 178);
    const auto *pfd_179 = buffer.data(pfd + 179);

#pragma omp simd aligned(t_255, t_256, t_257, t_258, t_259, pb_x, pc_x, pc_y, pdf0_176, \
                         pdf0_177, pdf0_179, pdd_107, pdf1_176, pdf1_177, pdf1_179, \
                         pfd_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_0 * pdd_107[k]
                   + f_4 * pc_x[k] * pfd_155[k];

        t_256[k] = pb_x[k] * pdf0_176[k]
                   - f_5 * pc_x[k] * pdf1_176[k];

        t_257[k] = pb_x[k] * pdf0_177[k]
                   - f_5 * pc_x[k] * pdf1_177[k];

        t_258[k] = f_4 * pc_y[k] * pfd_155[k];

        t_259[k] = pb_x[k] * pdf0_179[k]
                   - f_5 * pc_x[k] * pdf1_179[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pa_z, pc_x, pc_y, pc_z, sff0_60, sff1_60, \
                         pdd_90, pfp0_80, pfp1_80, pfd_156, pfd_158, \
                         pfd_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = pa_z[k] * sff0_60[k]
                   - f_5 * pc_z[k] * sff1_60[k];

        t_261[k] = f_1 * pdd_90[k]
                   + f_4 * pc_y[k] * pfd_156[k];

        t_262[k] = f_7 * pfp0_80[k]
                   - f_8 * pfp1_80[k]
                   + f_4 * pc_x[k] * pfd_158[k];

        t_263[k] = f_4 * pc_x[k] * pfd_159[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, pa_z, pc_x, pc_z, sff0_66, sff0_67, \
                         sfd_39, sff1_66, sff1_67, pfd_160, pfd_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_4 * pc_x[k] * pfd_160[k];

        t_265[k] = f_4 * pc_x[k] * pfd_161[k];

        t_266[k] = pa_z[k] * sff0_66[k]
                   - f_5 * pc_z[k] * sff1_66[k];

        t_267[k] = pa_z[k] * sff0_67[k]
                   + f_0 * sfd_39[k]
                   - f_5 * pc_z[k] * sff1_67[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, pa_z, pc_x, pc_y, pc_z, sff0_69, sfd_41, \
                         sff1_69, pdd_95, pfp0_81, pfp1_81, pfd_161, \
                         pfd_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_1 * pdd_95[k]
                   + f_4 * pc_y[k] * pfd_161[k];

        t_269[k] = pa_z[k] * sff0_69[k]
                   + f_1 * sfd_41[k]
                   - f_5 * pc_z[k] * sff1_69[k];

        t_270[k] = f_2 * pfp0_81[k]
                   - f_3 * pfp1_81[k]
                   + f_4 * pc_x[k] * pfd_162[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, t_275, pc_x, pc_y, pdd_96, pfp0_83, \
                         pfp1_83, pfd_162, pfd_164, pfd_165, pfd_166, \
                         pfd_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_6 * pdd_96[k]
                   + f_4 * pc_y[k] * pfd_162[k];

        t_272[k] = f_7 * pfp0_83[k]
                   - f_8 * pfp1_83[k]
                   + f_4 * pc_x[k] * pfd_164[k];

        t_273[k] = f_4 * pc_x[k] * pfd_165[k];

        t_274[k] = f_4 * pc_x[k] * pfd_166[k];

        t_275[k] = f_4 * pc_x[k] * pfd_167[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, pc_y, pdd_99, pdd_100, pdd_101, pfp0_82, \
                         pfp0_83, pfp1_82, pfp1_83, pfd_165, pfd_166, \
                         pfd_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_6 * pdd_99[k]
                   + f_2 * pfp0_82[k]
                   - f_3 * pfp1_82[k]
                   + f_4 * pc_y[k] * pfd_165[k];

        t_277[k] = f_6 * pdd_100[k]
                   + f_7 * pfp0_83[k]
                   - f_8 * pfp1_83[k]
                   + f_4 * pc_y[k] * pfd_166[k];

        t_278[k] = f_6 * pdd_101[k]
                   + f_4 * pc_y[k] * pfd_167[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, t_282, pb_y, pc_y, ppf0_89, ppf1_89, pdf0_169, \
                         pdf0_170, pdf0_172, pdd_102, pdf1_169, pdf1_170, pdf1_172, \
                         pfd_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = f_9 * ppf0_89[k]
                   - f_10 * ppf1_89[k]
                   + pb_y[k] * pdf0_169[k]
                   - f_5 * pc_y[k] * pdf1_169[k];

        t_280[k] = pb_y[k] * pdf0_170[k]
                   - f_5 * pc_y[k] * pdf1_170[k];

        t_281[k] = f_0 * pdd_102[k]
                   + f_4 * pc_y[k] * pfd_168[k];

        t_282[k] = pb_y[k] * pdf0_172[k]
                   - f_5 * pc_y[k] * pdf1_172[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, pb_y, pc_x, pc_y, pdf0_176, pdd_105, \
                         pdf1_176, pfd_171, pfd_172, pfd_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_4 * pc_x[k] * pfd_171[k];

        t_284[k] = f_4 * pc_x[k] * pfd_172[k];

        t_285[k] = f_4 * pc_x[k] * pfd_173[k];

        t_286[k] = pb_y[k] * pdf0_176[k]
                   + f_1 * pdd_105[k]
                   - f_5 * pc_y[k] * pdf1_176[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, pb_y, pc_y, pdf0_177, pdf0_179, pdd_106, \
                         pdd_107, pdf1_177, pdf1_179, pfd_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = pb_y[k] * pdf0_177[k]
                   + f_6 * pdd_106[k]
                   - f_5 * pc_y[k] * pdf1_177[k];

        t_288[k] = f_0 * pdd_107[k]
                   + f_4 * pc_y[k] * pfd_173[k];

        t_289[k] = pb_y[k] * pdf0_179[k]
                   - f_5 * pc_y[k] * pdf1_179[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, pc_x, pc_y, pfp0_87, pfp0_89, \
                         pfp1_87, pfp1_89, pfd_174, pfd_176, pfd_177, \
                         pfd_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = f_2 * pfp0_87[k]
                   - f_3 * pfp1_87[k]
                   + f_4 * pc_x[k] * pfd_174[k];

        t_291[k] = f_4 * pc_y[k] * pfd_174[k];

        t_292[k] = f_7 * pfp0_89[k]
                   - f_8 * pfp1_89[k]
                   + f_4 * pc_x[k] * pfd_176[k];

        t_293[k] = f_4 * pc_x[k] * pfd_177[k];

        t_294[k] = f_4 * pc_x[k] * pfd_178[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, pc_x, pc_y, pfp0_88, pfp0_89, pfp1_88, \
                         pfp1_89, pfd_177, pfd_178, pfd_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = f_4 * pc_x[k] * pfd_179[k];

        t_296[k] = f_2 * pfp0_88[k]
                   - f_3 * pfp1_88[k]
                   + f_4 * pc_y[k] * pfd_177[k];

        t_297[k] = f_7 * pfp0_89[k]
                   - f_8 * pfp1_89[k]
                   + f_4 * pc_y[k] * pfd_178[k];

        t_298[k] = f_4 * pc_y[k] * pfd_179[k];
    }

#pragma omp simd aligned(t_299, pc_z, sfd_59, pdd_107, pfp0_89, pfp1_89, \
                         pfd_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = f_0 * sfd_59[k]
                   + f_1 * pdd_107[k]
                   + f_2 * pfp0_89[k]
                   - f_3 * pfp1_89[k]
                   + f_4 * pc_z[k] * pfd_179[k];
    }
}

auto
compute_prim_pff_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t sff0,
                                                   const size_t sfd, const size_t sff1,
                                                   const size_t ppf0, const size_t ppf1,
                                                   const size_t pdf0, const size_t pdd,
                                                   const size_t pdf1, const size_t pfp0,
                                                   const size_t pfp1, const size_t pfd,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_pff_three_center_electron_repulsion_0_piece0(buffer, target, pa, pb, pc, sff0,
                                                              sfd, sff1, ppf0, ppf1, pdf0, pdd,
                                                              pdf1, pfp0, pfp1, pfd, ncols,
                                                              gamma, p, q);

    compute_prim_pff_three_center_electron_repulsion_0_piece1(buffer, target, pa, pb, pc, sff0,
                                                              sfd, sff1, ppf0, ppf1, pdf0, pdd,
                                                              pdf1, pfp0, pfp1, pfd, ncols,
                                                              gamma, p, q);

    compute_prim_pff_three_center_electron_repulsion_0_piece2(buffer, target, pa, pb, pc, sff0,
                                                              sfd, sff1, ppf0, ppf1, pdf0, pdd,
                                                              pdf1, pfp0, pfp1, pfd, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
