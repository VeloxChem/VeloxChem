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


#include "SimdThreeCenterElectronRepulsionVrrRecSID.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_sid_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shd0,
                                                          const size_t shp, const size_t shd1,
                                                          const size_t sis0, const size_t sis1,
                                                          const size_t sip, const size_t ncols,
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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *shd0_0 = buffer.data(shd0 + 0);
    const auto *shd0_3 = buffer.data(shd0 + 3);
    const auto *shd0_5 = buffer.data(shd0 + 5);
    const auto *shd0_9 = buffer.data(shd0 + 9);
    const auto *shd0_12 = buffer.data(shd0 + 12);
    const auto *shd0_17 = buffer.data(shd0 + 17);
    const auto *shd0_18 = buffer.data(shd0 + 18);
    const auto *shd0_21 = buffer.data(shd0 + 21);
    const auto *shd0_30 = buffer.data(shd0 + 30);
    const auto *shd0_35 = buffer.data(shd0 + 35);
    const auto *shd0_36 = buffer.data(shd0 + 36);
    const auto *shd0_39 = buffer.data(shd0 + 39);
    const auto *shd0_54 = buffer.data(shd0 + 54);
    const auto *shd0_59 = buffer.data(shd0 + 59);
    const auto *shd0_60 = buffer.data(shd0 + 60);
    const auto *shd0_84 = buffer.data(shd0 + 84);
    const auto *shd0_90 = buffer.data(shd0 + 90);
    const auto *shd0_93 = buffer.data(shd0 + 93);
    const auto *shd0_95 = buffer.data(shd0 + 95);
    const auto *shd0_99 = buffer.data(shd0 + 99);
    const auto *shd0_101 = buffer.data(shd0 + 101);
    const auto *shd0_102 = buffer.data(shd0 + 102);
    const auto *shd0_105 = buffer.data(shd0 + 105);
    const auto *shd0_107 = buffer.data(shd0 + 107);
    const auto *shd0_108 = buffer.data(shd0 + 108);
    const auto *shd0_111 = buffer.data(shd0 + 111);
    const auto *shd0_113 = buffer.data(shd0 + 113);
    const auto *shd0_117 = buffer.data(shd0 + 117);
    const auto *shd0_119 = buffer.data(shd0 + 119);
    const auto *shd0_120 = buffer.data(shd0 + 120);

    const auto *shp_0 = buffer.data(shp + 0);
    const auto *shp_1 = buffer.data(shp + 1);
    const auto *shp_2 = buffer.data(shp + 2);
    const auto *shp_4 = buffer.data(shp + 4);
    const auto *shp_5 = buffer.data(shp + 5);
    const auto *shp_7 = buffer.data(shp + 7);
    const auto *shp_8 = buffer.data(shp + 8);
    const auto *shp_9 = buffer.data(shp + 9);
    const auto *shp_10 = buffer.data(shp + 10);
    const auto *shp_11 = buffer.data(shp + 11);
    const auto *shp_13 = buffer.data(shp + 13);
    const auto *shp_14 = buffer.data(shp + 14);
    const auto *shp_15 = buffer.data(shp + 15);
    const auto *shp_16 = buffer.data(shp + 16);
    const auto *shp_17 = buffer.data(shp + 17);
    const auto *shp_18 = buffer.data(shp + 18);
    const auto *shp_19 = buffer.data(shp + 19);
    const auto *shp_20 = buffer.data(shp + 20);
    const auto *shp_22 = buffer.data(shp + 22);
    const auto *shp_23 = buffer.data(shp + 23);
    const auto *shp_25 = buffer.data(shp + 25);
    const auto *shp_26 = buffer.data(shp + 26);
    const auto *shp_27 = buffer.data(shp + 27);
    const auto *shp_28 = buffer.data(shp + 28);
    const auto *shp_29 = buffer.data(shp + 29);
    const auto *shp_30 = buffer.data(shp + 30);
    const auto *shp_31 = buffer.data(shp + 31);
    const auto *shp_32 = buffer.data(shp + 32);
    const auto *shp_34 = buffer.data(shp + 34);
    const auto *shp_35 = buffer.data(shp + 35);
    const auto *shp_36 = buffer.data(shp + 36);
    const auto *shp_37 = buffer.data(shp + 37);
    const auto *shp_38 = buffer.data(shp + 38);
    const auto *shp_40 = buffer.data(shp + 40);
    const auto *shp_41 = buffer.data(shp + 41);
    const auto *shp_42 = buffer.data(shp + 42);
    const auto *shp_43 = buffer.data(shp + 43);
    const auto *shp_44 = buffer.data(shp + 44);
    const auto *shp_45 = buffer.data(shp + 45);
    const auto *shp_46 = buffer.data(shp + 46);
    const auto *shp_47 = buffer.data(shp + 47);
    const auto *shp_49 = buffer.data(shp + 49);
    const auto *shp_50 = buffer.data(shp + 50);
    const auto *shp_51 = buffer.data(shp + 51);
    const auto *shp_52 = buffer.data(shp + 52);
    const auto *shp_53 = buffer.data(shp + 53);
    const auto *shp_54 = buffer.data(shp + 54);
    const auto *shp_55 = buffer.data(shp + 55);
    const auto *shp_56 = buffer.data(shp + 56);
    const auto *shp_58 = buffer.data(shp + 58);
    const auto *shp_59 = buffer.data(shp + 59);
    const auto *shp_60 = buffer.data(shp + 60);

    const auto *shd1_0 = buffer.data(shd1 + 0);
    const auto *shd1_3 = buffer.data(shd1 + 3);
    const auto *shd1_5 = buffer.data(shd1 + 5);
    const auto *shd1_9 = buffer.data(shd1 + 9);
    const auto *shd1_12 = buffer.data(shd1 + 12);
    const auto *shd1_17 = buffer.data(shd1 + 17);
    const auto *shd1_18 = buffer.data(shd1 + 18);
    const auto *shd1_21 = buffer.data(shd1 + 21);
    const auto *shd1_30 = buffer.data(shd1 + 30);
    const auto *shd1_35 = buffer.data(shd1 + 35);
    const auto *shd1_36 = buffer.data(shd1 + 36);
    const auto *shd1_39 = buffer.data(shd1 + 39);
    const auto *shd1_54 = buffer.data(shd1 + 54);
    const auto *shd1_59 = buffer.data(shd1 + 59);
    const auto *shd1_60 = buffer.data(shd1 + 60);
    const auto *shd1_84 = buffer.data(shd1 + 84);
    const auto *shd1_90 = buffer.data(shd1 + 90);
    const auto *shd1_93 = buffer.data(shd1 + 93);
    const auto *shd1_95 = buffer.data(shd1 + 95);
    const auto *shd1_99 = buffer.data(shd1 + 99);
    const auto *shd1_101 = buffer.data(shd1 + 101);
    const auto *shd1_102 = buffer.data(shd1 + 102);
    const auto *shd1_105 = buffer.data(shd1 + 105);
    const auto *shd1_107 = buffer.data(shd1 + 107);
    const auto *shd1_108 = buffer.data(shd1 + 108);
    const auto *shd1_111 = buffer.data(shd1 + 111);
    const auto *shd1_113 = buffer.data(shd1 + 113);
    const auto *shd1_117 = buffer.data(shd1 + 117);
    const auto *shd1_119 = buffer.data(shd1 + 119);
    const auto *shd1_120 = buffer.data(shd1 + 120);

    const auto *sis0_0 = buffer.data(sis0 + 0);
    const auto *sis0_1 = buffer.data(sis0 + 1);
    const auto *sis0_2 = buffer.data(sis0 + 2);
    const auto *sis0_3 = buffer.data(sis0 + 3);
    const auto *sis0_5 = buffer.data(sis0 + 5);
    const auto *sis0_6 = buffer.data(sis0 + 6);
    const auto *sis0_7 = buffer.data(sis0 + 7);
    const auto *sis0_8 = buffer.data(sis0 + 8);
    const auto *sis0_9 = buffer.data(sis0 + 9);
    const auto *sis0_10 = buffer.data(sis0 + 10);
    const auto *sis0_11 = buffer.data(sis0 + 11);
    const auto *sis0_12 = buffer.data(sis0 + 12);
    const auto *sis0_13 = buffer.data(sis0 + 13);
    const auto *sis0_14 = buffer.data(sis0 + 14);

    const auto *sis1_0 = buffer.data(sis1 + 0);
    const auto *sis1_1 = buffer.data(sis1 + 1);
    const auto *sis1_2 = buffer.data(sis1 + 2);
    const auto *sis1_3 = buffer.data(sis1 + 3);
    const auto *sis1_5 = buffer.data(sis1 + 5);
    const auto *sis1_6 = buffer.data(sis1 + 6);
    const auto *sis1_7 = buffer.data(sis1 + 7);
    const auto *sis1_8 = buffer.data(sis1 + 8);
    const auto *sis1_9 = buffer.data(sis1 + 9);
    const auto *sis1_10 = buffer.data(sis1 + 10);
    const auto *sis1_11 = buffer.data(sis1 + 11);
    const auto *sis1_12 = buffer.data(sis1 + 12);
    const auto *sis1_13 = buffer.data(sis1 + 13);
    const auto *sis1_14 = buffer.data(sis1 + 14);

    const auto *sip_0 = buffer.data(sip + 0);
    const auto *sip_1 = buffer.data(sip + 1);
    const auto *sip_2 = buffer.data(sip + 2);
    const auto *sip_4 = buffer.data(sip + 4);
    const auto *sip_5 = buffer.data(sip + 5);
    const auto *sip_7 = buffer.data(sip + 7);
    const auto *sip_8 = buffer.data(sip + 8);
    const auto *sip_9 = buffer.data(sip + 9);
    const auto *sip_10 = buffer.data(sip + 10);
    const auto *sip_11 = buffer.data(sip + 11);
    const auto *sip_13 = buffer.data(sip + 13);
    const auto *sip_14 = buffer.data(sip + 14);
    const auto *sip_15 = buffer.data(sip + 15);
    const auto *sip_16 = buffer.data(sip + 16);
    const auto *sip_17 = buffer.data(sip + 17);
    const auto *sip_18 = buffer.data(sip + 18);
    const auto *sip_19 = buffer.data(sip + 19);
    const auto *sip_20 = buffer.data(sip + 20);
    const auto *sip_22 = buffer.data(sip + 22);
    const auto *sip_23 = buffer.data(sip + 23);
    const auto *sip_25 = buffer.data(sip + 25);
    const auto *sip_26 = buffer.data(sip + 26);
    const auto *sip_27 = buffer.data(sip + 27);
    const auto *sip_28 = buffer.data(sip + 28);
    const auto *sip_29 = buffer.data(sip + 29);
    const auto *sip_30 = buffer.data(sip + 30);
    const auto *sip_31 = buffer.data(sip + 31);
    const auto *sip_32 = buffer.data(sip + 32);
    const auto *sip_34 = buffer.data(sip + 34);
    const auto *sip_35 = buffer.data(sip + 35);
    const auto *sip_36 = buffer.data(sip + 36);
    const auto *sip_37 = buffer.data(sip + 37);
    const auto *sip_38 = buffer.data(sip + 38);
    const auto *sip_40 = buffer.data(sip + 40);
    const auto *sip_41 = buffer.data(sip + 41);
    const auto *sip_42 = buffer.data(sip + 42);
    const auto *sip_43 = buffer.data(sip + 43);
    const auto *sip_44 = buffer.data(sip + 44);
    const auto *sip_46 = buffer.data(sip + 46);
    const auto *sip_47 = buffer.data(sip + 47);
    const auto *sip_49 = buffer.data(sip + 49);
    const auto *sip_50 = buffer.data(sip + 50);
    const auto *sip_52 = buffer.data(sip + 52);
    const auto *sip_53 = buffer.data(sip + 53);
    const auto *sip_55 = buffer.data(sip + 55);
    const auto *sip_56 = buffer.data(sip + 56);
    const auto *sip_58 = buffer.data(sip + 58);
    const auto *sip_59 = buffer.data(sip + 59);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, shp_0, shp_1, shp_2, sis0_0, \
                         sis1_0, sip_0, sip_1, sip_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * shp_0[k]
                 + f_1 * sis0_0[k]
                 - f_2 * sis1_0[k]
                 + f_3 * pc_x[k] * sip_0[k];

        t_1[k] = f_0 * shp_1[k]
                 + f_3 * pc_x[k] * sip_1[k];

        t_2[k] = f_0 * shp_2[k]
                 + f_3 * pc_x[k] * sip_2[k];

        t_3[k] = f_1 * sis0_0[k]
                 - f_2 * sis1_0[k]
                 + f_3 * pc_y[k] * sip_1[k];

        t_4[k] = f_3 * pc_y[k] * sip_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pb_y, pc_x, pc_y, pc_z, shd0_0, shp_4, shd1_0, sis0_0, \
                         sis1_0, sip_2, sip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * sis0_0[k]
                 - f_2 * sis1_0[k]
                 + f_3 * pc_z[k] * sip_2[k];

        t_6[k] = pb_y[k] * shd0_0[k]
                 - f_4 * pc_y[k] * shd1_0[k];

        t_7[k] = f_5 * shp_4[k]
                 + f_3 * pc_x[k] * sip_4[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_y, pc_x, pc_y, shd0_5, shp_1, shp_2, shp_5, \
                         shd1_5, sis0_1, sis1_1, sip_4, sip_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_5 * shp_5[k]
                 + f_3 * pc_x[k] * sip_5[k];

        t_9[k] = f_6 * shp_1[k]
                 + f_1 * sis0_1[k]
                 - f_2 * sis1_1[k]
                 + f_3 * pc_y[k] * sip_4[k];

        t_10[k] = f_6 * shp_2[k]
                  + f_3 * pc_y[k] * sip_5[k];

        t_11[k] = pb_y[k] * shd0_5[k]
                  - f_4 * pc_y[k] * shd1_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pb_z, pc_x, pc_z, shd0_0, shd0_3, shp_7, \
                         shp_8, shd1_0, shd1_3, sip_7, sip_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pb_z[k] * shd0_0[k]
                  - f_4 * pc_z[k] * shd1_0[k];

        t_13[k] = f_5 * shp_7[k]
                  + f_3 * pc_x[k] * sip_7[k];

        t_14[k] = f_5 * shp_8[k]
                  + f_3 * pc_x[k] * sip_8[k];

        t_15[k] = pb_z[k] * shd0_3[k]
                  - f_4 * pc_z[k] * shd1_3[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pc_x, pc_y, pc_z, shp_2, shp_9, sis0_2, sis0_3, \
                         sis1_2, sis1_3, sip_8, sip_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_y[k] * sip_8[k];

        t_17[k] = f_6 * shp_2[k]
                  + f_1 * sis0_2[k]
                  - f_2 * sis1_2[k]
                  + f_3 * pc_z[k] * sip_8[k];

        t_18[k] = f_7 * shp_9[k]
                  + f_1 * sis0_3[k]
                  - f_2 * sis1_3[k]
                  + f_3 * pc_x[k] * sip_9[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pc_x, pc_y, pc_z, shp_4, shp_5, shp_10, \
                         shp_11, sis0_3, sis1_3, sip_10, sip_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_7 * shp_10[k]
                  + f_3 * pc_x[k] * sip_10[k];

        t_20[k] = f_7 * shp_11[k]
                  + f_3 * pc_x[k] * sip_11[k];

        t_21[k] = f_8 * shp_4[k]
                  + f_1 * sis0_3[k]
                  - f_2 * sis1_3[k]
                  + f_3 * pc_y[k] * sip_10[k];

        t_22[k] = f_8 * shp_5[k]
                  + f_3 * pc_y[k] * sip_11[k];

        t_23[k] = f_1 * sis0_3[k]
                  - f_2 * sis1_3[k]
                  + f_3 * pc_z[k] * sip_11[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pb_y, pc_x, pc_y, shd0_12, shp_13, shp_14, shd1_12, \
                         sip_13, sip_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pb_y[k] * shd0_12[k]
                  - f_4 * pc_y[k] * shd1_12[k];

        t_25[k] = f_7 * shp_13[k]
                  + f_3 * pc_x[k] * sip_13[k];

        t_26[k] = f_7 * shp_14[k]
                  + f_3 * pc_x[k] * sip_14[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pb_y, pb_z, pc_y, pc_z, shd0_9, shd0_17, shp_8, \
                         shd1_9, shd1_17, sip_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pb_z[k] * shd0_9[k]
                  - f_4 * pc_z[k] * shd1_9[k];

        t_28[k] = f_6 * shp_8[k]
                  + f_3 * pc_y[k] * sip_14[k];

        t_29[k] = pb_y[k] * shd0_17[k]
                  - f_4 * pc_y[k] * shd1_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pc_x, pc_y, shp_15, shp_16, shp_17, \
                         sis0_5, sis1_5, sip_15, sip_16, sip_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_7 * shp_15[k]
                  + f_1 * sis0_5[k]
                  - f_2 * sis1_5[k]
                  + f_3 * pc_x[k] * sip_15[k];

        t_31[k] = f_7 * shp_16[k]
                  + f_3 * pc_x[k] * sip_16[k];

        t_32[k] = f_7 * shp_17[k]
                  + f_3 * pc_x[k] * sip_17[k];

        t_33[k] = f_1 * sis0_5[k]
                  - f_2 * sis1_5[k]
                  + f_3 * pc_y[k] * sip_16[k];

        t_34[k] = f_3 * pc_y[k] * sip_17[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, pc_x, pc_z, shp_8, shp_18, shp_19, sis0_5, sis0_6, \
                         sis1_5, sis1_6, sip_17, sip_18, sip_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_8 * shp_8[k]
                  + f_1 * sis0_5[k]
                  - f_2 * sis1_5[k]
                  + f_3 * pc_z[k] * sip_17[k];

        t_36[k] = f_9 * shp_18[k]
                  + f_1 * sis0_6[k]
                  - f_2 * sis1_6[k]
                  + f_3 * pc_x[k] * sip_18[k];

        t_37[k] = f_9 * shp_19[k]
                  + f_3 * pc_x[k] * sip_19[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pc_x, pc_y, pc_z, shp_10, shp_11, shp_20, \
                         sis0_6, sis1_6, sip_19, sip_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = f_9 * shp_20[k]
                  + f_3 * pc_x[k] * sip_20[k];

        t_39[k] = f_9 * shp_10[k]
                  + f_1 * sis0_6[k]
                  - f_2 * sis1_6[k]
                  + f_3 * pc_y[k] * sip_19[k];

        t_40[k] = f_9 * shp_11[k]
                  + f_3 * pc_y[k] * sip_20[k];

        t_41[k] = f_1 * sis0_6[k]
                  - f_2 * sis1_6[k]
                  + f_3 * pc_z[k] * sip_20[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pb_z, pc_x, pc_z, shd0_18, shd0_21, shp_22, \
                         shp_23, shd1_18, shd1_21, sip_22, sip_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_z[k] * shd0_18[k]
                  - f_4 * pc_z[k] * shd1_18[k];

        t_43[k] = f_9 * shp_22[k]
                  + f_3 * pc_x[k] * sip_22[k];

        t_44[k] = f_9 * shp_23[k]
                  + f_3 * pc_x[k] * sip_23[k];

        t_45[k] = pb_z[k] * shd0_21[k]
                  - f_4 * pc_z[k] * shd1_21[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pb_y, pc_y, pc_z, shd0_30, shp_11, shp_14, shd1_30, \
                         sis0_7, sis1_7, sip_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_8 * shp_14[k]
                  + f_3 * pc_y[k] * sip_23[k];

        t_47[k] = f_6 * shp_11[k]
                  + f_1 * sis0_7[k]
                  - f_2 * sis1_7[k]
                  + f_3 * pc_z[k] * sip_23[k];

        t_48[k] = pb_y[k] * shd0_30[k]
                  - f_4 * pc_y[k] * shd1_30[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pc_x, pc_y, shp_16, shp_17, shp_25, shp_26, \
                         sis0_8, sis1_8, sip_25, sip_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_9 * shp_25[k]
                  + f_3 * pc_x[k] * sip_25[k];

        t_50[k] = f_9 * shp_26[k]
                  + f_3 * pc_x[k] * sip_26[k];

        t_51[k] = f_6 * shp_16[k]
                  + f_1 * sis0_8[k]
                  - f_2 * sis1_8[k]
                  + f_3 * pc_y[k] * sip_25[k];

        t_52[k] = f_6 * shp_17[k]
                  + f_3 * pc_y[k] * sip_26[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_y, pc_x, pc_y, shd0_35, shp_27, shp_28, shd1_35, \
                         sis0_9, sis1_9, sip_27, sip_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = pb_y[k] * shd0_35[k]
                  - f_4 * pc_y[k] * shd1_35[k];

        t_54[k] = f_9 * shp_27[k]
                  + f_1 * sis0_9[k]
                  - f_2 * sis1_9[k]
                  + f_3 * pc_x[k] * sip_27[k];

        t_55[k] = f_9 * shp_28[k]
                  + f_3 * pc_x[k] * sip_28[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, pc_x, pc_y, pc_z, shp_17, shp_29, sis0_9, \
                         sis1_9, sip_28, sip_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_9 * shp_29[k]
                  + f_3 * pc_x[k] * sip_29[k];

        t_57[k] = f_1 * sis0_9[k]
                  - f_2 * sis1_9[k]
                  + f_3 * pc_y[k] * sip_28[k];

        t_58[k] = f_3 * pc_y[k] * sip_29[k];

        t_59[k] = f_9 * shp_17[k]
                  + f_1 * sis0_9[k]
                  - f_2 * sis1_9[k]
                  + f_3 * pc_z[k] * sip_29[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, t_63, pc_x, pc_y, shp_19, shp_30, shp_31, shp_32, \
                         sis0_10, sis1_10, sip_30, sip_31, sip_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_8 * shp_30[k]
                  + f_1 * sis0_10[k]
                  - f_2 * sis1_10[k]
                  + f_3 * pc_x[k] * sip_30[k];

        t_61[k] = f_8 * shp_31[k]
                  + f_3 * pc_x[k] * sip_31[k];

        t_62[k] = f_8 * shp_32[k]
                  + f_3 * pc_x[k] * sip_32[k];

        t_63[k] = f_7 * shp_19[k]
                  + f_1 * sis0_10[k]
                  - f_2 * sis1_10[k]
                  + f_3 * pc_y[k] * sip_31[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, t_67, pb_z, pc_x, pc_y, pc_z, shd0_36, shp_20, \
                         shp_34, shd1_36, sis0_10, sis1_10, sip_32, \
                         sip_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_7 * shp_20[k]
                  + f_3 * pc_y[k] * sip_32[k];

        t_65[k] = f_1 * sis0_10[k]
                  - f_2 * sis1_10[k]
                  + f_3 * pc_z[k] * sip_32[k];

        t_66[k] = pb_z[k] * shd0_36[k]
                  - f_4 * pc_z[k] * shd1_36[k];

        t_67[k] = f_8 * shp_34[k]
                  + f_3 * pc_x[k] * sip_34[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, pb_z, pc_x, pc_y, pc_z, shd0_39, shp_20, \
                         shp_23, shp_35, shd1_39, sis0_11, sis1_11, \
                         sip_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_8 * shp_35[k]
                  + f_3 * pc_x[k] * sip_35[k];

        t_69[k] = pb_z[k] * shd0_39[k]
                  - f_4 * pc_z[k] * shd1_39[k];

        t_70[k] = f_9 * shp_23[k]
                  + f_3 * pc_y[k] * sip_35[k];

        t_71[k] = f_6 * shp_20[k]
                  + f_1 * sis0_11[k]
                  - f_2 * sis1_11[k]
                  + f_3 * pc_z[k] * sip_35[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, pc_x, pc_y, shp_25, shp_36, shp_37, shp_38, \
                         sis0_12, sis1_12, sip_36, sip_37, sip_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_8 * shp_36[k]
                  + f_1 * sis0_12[k]
                  - f_2 * sis1_12[k]
                  + f_3 * pc_x[k] * sip_36[k];

        t_73[k] = f_8 * shp_37[k]
                  + f_3 * pc_x[k] * sip_37[k];

        t_74[k] = f_8 * shp_38[k]
                  + f_3 * pc_x[k] * sip_38[k];

        t_75[k] = f_8 * shp_25[k]
                  + f_1 * sis0_12[k]
                  - f_2 * sis1_12[k]
                  + f_3 * pc_y[k] * sip_37[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, pb_y, pc_y, pc_z, shd0_54, shp_23, shp_26, shd1_54, \
                         sis0_12, sis1_12, sip_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_8 * shp_26[k]
                  + f_3 * pc_y[k] * sip_38[k];

        t_77[k] = f_8 * shp_23[k]
                  + f_1 * sis0_12[k]
                  - f_2 * sis1_12[k]
                  + f_3 * pc_z[k] * sip_38[k];

        t_78[k] = pb_y[k] * shd0_54[k]
                  - f_4 * pc_y[k] * shd1_54[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pc_x, pc_y, shp_28, shp_29, shp_40, shp_41, \
                         sis0_13, sis1_13, sip_40, sip_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_8 * shp_40[k]
                  + f_3 * pc_x[k] * sip_40[k];

        t_80[k] = f_8 * shp_41[k]
                  + f_3 * pc_x[k] * sip_41[k];

        t_81[k] = f_6 * shp_28[k]
                  + f_1 * sis0_13[k]
                  - f_2 * sis1_13[k]
                  + f_3 * pc_y[k] * sip_40[k];

        t_82[k] = f_6 * shp_29[k]
                  + f_3 * pc_y[k] * sip_41[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, pb_y, pc_x, pc_y, shd0_59, shp_42, shp_43, shd1_59, \
                         sis0_14, sis1_14, sip_42, sip_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = pb_y[k] * shd0_59[k]
                  - f_4 * pc_y[k] * shd1_59[k];

        t_84[k] = f_8 * shp_42[k]
                  + f_1 * sis0_14[k]
                  - f_2 * sis1_14[k]
                  + f_3 * pc_x[k] * sip_42[k];

        t_85[k] = f_8 * shp_43[k]
                  + f_3 * pc_x[k] * sip_43[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pc_x, pc_y, pc_z, shp_29, shp_44, sis0_14, \
                         sis1_14, sip_43, sip_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_8 * shp_44[k]
                  + f_3 * pc_x[k] * sip_44[k];

        t_87[k] = f_1 * sis0_14[k]
                  - f_2 * sis1_14[k]
                  + f_3 * pc_y[k] * sip_43[k];

        t_88[k] = f_3 * pc_y[k] * sip_44[k];

        t_89[k] = f_7 * shp_29[k]
                  + f_1 * sis0_14[k]
                  - f_2 * sis1_14[k]
                  + f_3 * pc_z[k] * sip_44[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pb_x, pc_x, shd0_90, shd0_93, shp_45, shp_46, \
                         shp_47, shd1_90, shd1_93, sip_46, sip_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = pb_x[k] * shd0_90[k]
                  + f_8 * shp_45[k]
                  - f_4 * pc_x[k] * shd1_90[k];

        t_91[k] = f_6 * shp_46[k]
                  + f_3 * pc_x[k] * sip_46[k];

        t_92[k] = f_6 * shp_47[k]
                  + f_3 * pc_x[k] * sip_47[k];

        t_93[k] = pb_x[k] * shd0_93[k]
                  - f_4 * pc_x[k] * shd1_93[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pb_x, pb_z, pc_x, pc_y, pc_z, shd0_60, shd0_95, \
                         shp_32, shd1_60, shd1_95, sip_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_5 * shp_32[k]
                  + f_3 * pc_y[k] * sip_47[k];

        t_95[k] = pb_x[k] * shd0_95[k]
                  - f_4 * pc_x[k] * shd1_95[k];

        t_96[k] = pb_z[k] * shd0_60[k]
                  - f_4 * pc_z[k] * shd1_60[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pb_x, pc_x, pc_y, shd0_99, shp_35, shp_49, \
                         shp_50, shd1_99, sip_49, sip_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_6 * shp_49[k]
                  + f_3 * pc_x[k] * sip_49[k];

        t_98[k] = f_6 * shp_50[k]
                  + f_3 * pc_x[k] * sip_50[k];

        t_99[k] = pb_x[k] * shd0_99[k]
                  - f_4 * pc_x[k] * shd1_99[k];

        t_100[k] = f_7 * shp_35[k]
                   + f_3 * pc_y[k] * sip_50[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pb_x, pc_x, shd0_101, shd0_102, shp_51, \
                         shp_52, shp_53, shd1_101, shd1_102, sip_52, \
                         sip_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = pb_x[k] * shd0_101[k]
                   - f_4 * pc_x[k] * shd1_101[k];

        t_102[k] = pb_x[k] * shd0_102[k]
                   + f_8 * shp_51[k]
                   - f_4 * pc_x[k] * shd1_102[k];

        t_103[k] = f_6 * shp_52[k]
                   + f_3 * pc_x[k] * sip_52[k];

        t_104[k] = f_6 * shp_53[k]
                   + f_3 * pc_x[k] * sip_53[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pb_x, pc_x, pc_y, shd0_105, shd0_107, \
                         shd0_108, shp_38, shp_54, shd1_105, shd1_107, shd1_108, \
                         sip_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = pb_x[k] * shd0_105[k]
                   - f_4 * pc_x[k] * shd1_105[k];

        t_106[k] = f_9 * shp_38[k]
                   + f_3 * pc_y[k] * sip_53[k];

        t_107[k] = pb_x[k] * shd0_107[k]
                   - f_4 * pc_x[k] * shd1_107[k];

        t_108[k] = pb_x[k] * shd0_108[k]
                   + f_8 * shp_54[k]
                   - f_4 * pc_x[k] * shd1_108[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pb_x, pc_x, pc_y, shd0_111, shp_41, \
                         shp_55, shp_56, shd1_111, sip_55, sip_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_6 * shp_55[k]
                   + f_3 * pc_x[k] * sip_55[k];

        t_110[k] = f_6 * shp_56[k]
                   + f_3 * pc_x[k] * sip_56[k];

        t_111[k] = pb_x[k] * shd0_111[k]
                   - f_4 * pc_x[k] * shd1_111[k];

        t_112[k] = f_8 * shp_41[k]
                   + f_3 * pc_y[k] * sip_56[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pb_x, pb_y, pc_x, pc_y, shd0_84, \
                         shd0_113, shp_58, shp_59, shd1_84, shd1_113, sip_58, \
                         sip_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = pb_x[k] * shd0_113[k]
                   - f_4 * pc_x[k] * shd1_113[k];

        t_114[k] = pb_y[k] * shd0_84[k]
                   - f_4 * pc_y[k] * shd1_84[k];

        t_115[k] = f_6 * shp_58[k]
                   + f_3 * pc_x[k] * sip_58[k];

        t_116[k] = f_6 * shp_59[k]
                   + f_3 * pc_x[k] * sip_59[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pb_x, pc_x, pc_y, shd0_117, shd0_119, \
                         shd0_120, shp_44, shp_60, shd1_117, shd1_119, shd1_120, \
                         sip_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = pb_x[k] * shd0_117[k]
                   - f_4 * pc_x[k] * shd1_117[k];

        t_118[k] = f_6 * shp_44[k]
                   + f_3 * pc_y[k] * sip_59[k];

        t_119[k] = pb_x[k] * shd0_119[k]
                   - f_4 * pc_x[k] * shd1_119[k];

        t_120[k] = pb_x[k] * shd0_120[k]
                   + f_8 * shp_60[k]
                   - f_4 * pc_x[k] * shd1_120[k];
    }
}

static auto
compute_prim_sid_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t shd0,
                                                          const size_t shp, const size_t shd1,
                                                          const size_t sis0, const size_t sis1,
                                                          const size_t sip, const size_t ncols,
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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *shd0_90 = buffer.data(shd0 + 90);
    const auto *shd0_93 = buffer.data(shd0 + 93);
    const auto *shd0_120 = buffer.data(shd0 + 120);
    const auto *shd0_123 = buffer.data(shd0 + 123);
    const auto *shd0_125 = buffer.data(shd0 + 125);

    const auto *shp_46 = buffer.data(shp + 46);
    const auto *shp_47 = buffer.data(shp + 47);
    const auto *shp_50 = buffer.data(shp + 50);
    const auto *shp_52 = buffer.data(shp + 52);
    const auto *shp_53 = buffer.data(shp + 53);
    const auto *shp_55 = buffer.data(shp + 55);
    const auto *shp_56 = buffer.data(shp + 56);
    const auto *shp_58 = buffer.data(shp + 58);
    const auto *shp_59 = buffer.data(shp + 59);
    const auto *shp_61 = buffer.data(shp + 61);
    const auto *shp_62 = buffer.data(shp + 62);

    const auto *shd1_90 = buffer.data(shd1 + 90);
    const auto *shd1_93 = buffer.data(shd1 + 93);
    const auto *shd1_120 = buffer.data(shd1 + 120);
    const auto *shd1_123 = buffer.data(shd1 + 123);
    const auto *shd1_125 = buffer.data(shd1 + 125);

    const auto *sis0_21 = buffer.data(sis0 + 21);
    const auto *sis0_22 = buffer.data(sis0 + 22);
    const auto *sis0_23 = buffer.data(sis0 + 23);
    const auto *sis0_24 = buffer.data(sis0 + 24);
    const auto *sis0_25 = buffer.data(sis0 + 25);
    const auto *sis0_27 = buffer.data(sis0 + 27);

    const auto *sis1_21 = buffer.data(sis1 + 21);
    const auto *sis1_22 = buffer.data(sis1 + 22);
    const auto *sis1_23 = buffer.data(sis1 + 23);
    const auto *sis1_24 = buffer.data(sis1 + 24);
    const auto *sis1_25 = buffer.data(sis1 + 25);
    const auto *sis1_27 = buffer.data(sis1 + 27);

    const auto *sip_61 = buffer.data(sip + 61);
    const auto *sip_62 = buffer.data(sip + 62);
    const auto *sip_63 = buffer.data(sip + 63);
    const auto *sip_64 = buffer.data(sip + 64);
    const auto *sip_65 = buffer.data(sip + 65);
    const auto *sip_67 = buffer.data(sip + 67);
    const auto *sip_68 = buffer.data(sip + 68);
    const auto *sip_69 = buffer.data(sip + 69);
    const auto *sip_70 = buffer.data(sip + 70);
    const auto *sip_71 = buffer.data(sip + 71);
    const auto *sip_72 = buffer.data(sip + 72);
    const auto *sip_73 = buffer.data(sip + 73);
    const auto *sip_74 = buffer.data(sip + 74);
    const auto *sip_75 = buffer.data(sip + 75);
    const auto *sip_76 = buffer.data(sip + 76);
    const auto *sip_77 = buffer.data(sip + 77);
    const auto *sip_79 = buffer.data(sip + 79);
    const auto *sip_80 = buffer.data(sip + 80);
    const auto *sip_81 = buffer.data(sip + 81);
    const auto *sip_82 = buffer.data(sip + 82);
    const auto *sip_83 = buffer.data(sip + 83);

#pragma omp simd aligned(t_121, t_122, t_123, t_124, t_125, pb_x, pc_x, pc_y, shd0_123, \
                         shd0_125, shp_61, shp_62, shd1_123, shd1_125, sip_61, \
                         sip_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_6 * shp_61[k]
                   + f_3 * pc_x[k] * sip_61[k];

        t_122[k] = f_6 * shp_62[k]
                   + f_3 * pc_x[k] * sip_62[k];

        t_123[k] = pb_x[k] * shd0_123[k]
                   - f_4 * pc_x[k] * shd1_123[k];

        t_124[k] = f_3 * pc_y[k] * sip_62[k];

        t_125[k] = pb_x[k] * shd0_125[k]
                   - f_4 * pc_x[k] * shd1_125[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, t_130, t_131, pc_x, pc_y, pc_z, shp_46, \
                         shp_47, sis0_21, sis1_21, sip_63, sip_64, \
                         sip_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = f_1 * sis0_21[k]
                   - f_2 * sis1_21[k]
                   + f_3 * pc_x[k] * sip_63[k];

        t_127[k] = f_3 * pc_x[k] * sip_64[k];

        t_128[k] = f_3 * pc_x[k] * sip_65[k];

        t_129[k] = f_0 * shp_46[k]
                   + f_1 * sis0_21[k]
                   - f_2 * sis1_21[k]
                   + f_3 * pc_y[k] * sip_64[k];

        t_130[k] = f_0 * shp_47[k]
                   + f_3 * pc_y[k] * sip_65[k];

        t_131[k] = f_1 * sis0_21[k]
                   - f_2 * sis1_21[k]
                   + f_3 * pc_z[k] * sip_65[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, t_136, pb_z, pc_x, pc_y, pc_z, shd0_90, \
                         shd0_93, shp_50, shd1_90, shd1_93, sip_67, \
                         sip_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = pb_z[k] * shd0_90[k]
                   - f_4 * pc_z[k] * shd1_90[k];

        t_133[k] = f_3 * pc_x[k] * sip_67[k];

        t_134[k] = f_3 * pc_x[k] * sip_68[k];

        t_135[k] = pb_z[k] * shd0_93[k]
                   - f_4 * pc_z[k] * shd1_93[k];

        t_136[k] = f_5 * shp_50[k]
                   + f_3 * pc_y[k] * sip_68[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, t_140, pc_x, pc_z, shp_47, sis0_22, sis0_23, \
                         sis1_22, sis1_23, sip_68, sip_69, sip_70, \
                         sip_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_6 * shp_47[k]
                   + f_1 * sis0_22[k]
                   - f_2 * sis1_22[k]
                   + f_3 * pc_z[k] * sip_68[k];

        t_138[k] = f_1 * sis0_23[k]
                   - f_2 * sis1_23[k]
                   + f_3 * pc_x[k] * sip_69[k];

        t_139[k] = f_3 * pc_x[k] * sip_70[k];

        t_140[k] = f_3 * pc_x[k] * sip_71[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, pc_y, pc_z, shp_50, shp_52, shp_53, sis0_23, \
                         sis1_23, sip_70, sip_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_7 * shp_52[k]
                   + f_1 * sis0_23[k]
                   - f_2 * sis1_23[k]
                   + f_3 * pc_y[k] * sip_70[k];

        t_142[k] = f_7 * shp_53[k]
                   + f_3 * pc_y[k] * sip_71[k];

        t_143[k] = f_8 * shp_50[k]
                   + f_1 * sis0_23[k]
                   - f_2 * sis1_23[k]
                   + f_3 * pc_z[k] * sip_71[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, t_148, pc_x, pc_y, shp_55, shp_56, \
                         sis0_24, sis1_24, sip_72, sip_73, sip_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_1 * sis0_24[k]
                   - f_2 * sis1_24[k]
                   + f_3 * pc_x[k] * sip_72[k];

        t_145[k] = f_3 * pc_x[k] * sip_73[k];

        t_146[k] = f_3 * pc_x[k] * sip_74[k];

        t_147[k] = f_9 * shp_55[k]
                   + f_1 * sis0_24[k]
                   - f_2 * sis1_24[k]
                   + f_3 * pc_y[k] * sip_73[k];

        t_148[k] = f_9 * shp_56[k]
                   + f_3 * pc_y[k] * sip_74[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, t_152, pc_x, pc_z, shp_53, sis0_24, sis0_25, \
                         sis1_24, sis1_25, sip_74, sip_75, sip_76, \
                         sip_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_9 * shp_53[k]
                   + f_1 * sis0_24[k]
                   - f_2 * sis1_24[k]
                   + f_3 * pc_z[k] * sip_74[k];

        t_150[k] = f_1 * sis0_25[k]
                   - f_2 * sis1_25[k]
                   + f_3 * pc_x[k] * sip_75[k];

        t_151[k] = f_3 * pc_x[k] * sip_76[k];

        t_152[k] = f_3 * pc_x[k] * sip_77[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, pb_y, pc_y, pc_z, shd0_120, shp_56, \
                         shp_58, shp_59, shd1_120, sis0_25, sis1_25, sip_76, \
                         sip_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_8 * shp_58[k]
                   + f_1 * sis0_25[k]
                   - f_2 * sis1_25[k]
                   + f_3 * pc_y[k] * sip_76[k];

        t_154[k] = f_8 * shp_59[k]
                   + f_3 * pc_y[k] * sip_77[k];

        t_155[k] = f_7 * shp_56[k]
                   + f_1 * sis0_25[k]
                   - f_2 * sis1_25[k]
                   + f_3 * pc_z[k] * sip_77[k];

        t_156[k] = pb_y[k] * shd0_120[k]
                   - f_4 * pc_y[k] * shd1_120[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, pb_y, pc_x, pc_y, shd0_123, \
                         shd0_125, shp_61, shp_62, shd1_123, shd1_125, sip_79, \
                         sip_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_3 * pc_x[k] * sip_79[k];

        t_158[k] = f_3 * pc_x[k] * sip_80[k];

        t_159[k] = pb_y[k] * shd0_123[k]
                   + f_8 * shp_61[k]
                   - f_4 * pc_y[k] * shd1_123[k];

        t_160[k] = f_6 * shp_62[k]
                   + f_3 * pc_y[k] * sip_80[k];

        t_161[k] = pb_y[k] * shd0_125[k]
                   - f_4 * pc_y[k] * shd1_125[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, t_166, t_167, pc_x, pc_y, pc_z, shp_62, \
                         sis0_27, sis1_27, sip_81, sip_82, sip_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_1 * sis0_27[k]
                   - f_2 * sis1_27[k]
                   + f_3 * pc_x[k] * sip_81[k];

        t_163[k] = f_3 * pc_x[k] * sip_82[k];

        t_164[k] = f_3 * pc_x[k] * sip_83[k];

        t_165[k] = f_1 * sis0_27[k]
                   - f_2 * sis1_27[k]
                   + f_3 * pc_y[k] * sip_82[k];

        t_166[k] = f_3 * pc_y[k] * sip_83[k];

        t_167[k] = f_0 * shp_62[k]
                   + f_1 * sis0_27[k]
                   - f_2 * sis1_27[k]
                   + f_3 * pc_z[k] * sip_83[k];
    }
}

auto
compute_prim_sid_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t shd0, const size_t shp,
                                                   const size_t shd1, const size_t sis0,
                                                   const size_t sis1, const size_t sip,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_sid_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, shd0, shp,
                                                              shd1, sis0, sis1, sip, ncols,
                                                              gamma, p, q);

    compute_prim_sid_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, shd0, shp,
                                                              shd1, sis0, sis1, sip, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
