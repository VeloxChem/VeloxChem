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


#include "SimdThreeCenterElectronRepulsionVrrRecSKG.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_skg_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sig0,
                                                          const size_t sif, const size_t sig1,
                                                          const size_t skd0, const size_t skd1,
                                                          const size_t skf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.0 / q;
    const auto f_10 = 2.5 / q;
    const auto f_11 = 2.0 / q;
    const auto f_12 = 1.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sig0_0 = buffer.data(sig0 + 0);
    const auto *sig0_3 = buffer.data(sig0 + 3);
    const auto *sig0_5 = buffer.data(sig0 + 5);
    const auto *sig0_10 = buffer.data(sig0 + 10);
    const auto *sig0_14 = buffer.data(sig0 + 14);
    const auto *sig0_18 = buffer.data(sig0 + 18);
    const auto *sig0_25 = buffer.data(sig0 + 25);
    const auto *sig0_30 = buffer.data(sig0 + 30);
    const auto *sig0_35 = buffer.data(sig0 + 35);
    const auto *sig0_44 = buffer.data(sig0 + 44);
    const auto *sig0_45 = buffer.data(sig0 + 45);
    const auto *sig0_48 = buffer.data(sig0 + 48);
    const auto *sig0_55 = buffer.data(sig0 + 55);
    const auto *sig0_75 = buffer.data(sig0 + 75);
    const auto *sig0_78 = buffer.data(sig0 + 78);

    const auto *sif_0 = buffer.data(sif + 0);
    const auto *sif_1 = buffer.data(sif + 1);
    const auto *sif_2 = buffer.data(sif + 2);
    const auto *sif_3 = buffer.data(sif + 3);
    const auto *sif_5 = buffer.data(sif + 5);
    const auto *sif_6 = buffer.data(sif + 6);
    const auto *sif_7 = buffer.data(sif + 7);
    const auto *sif_8 = buffer.data(sif + 8);
    const auto *sif_9 = buffer.data(sif + 9);
    const auto *sif_10 = buffer.data(sif + 10);
    const auto *sif_12 = buffer.data(sif + 12);
    const auto *sif_16 = buffer.data(sif + 16);
    const auto *sif_17 = buffer.data(sif + 17);
    const auto *sif_18 = buffer.data(sif + 18);
    const auto *sif_19 = buffer.data(sif + 19);
    const auto *sif_20 = buffer.data(sif + 20);
    const auto *sif_22 = buffer.data(sif + 22);
    const auto *sif_26 = buffer.data(sif + 26);
    const auto *sif_27 = buffer.data(sif + 27);
    const auto *sif_28 = buffer.data(sif + 28);
    const auto *sif_29 = buffer.data(sif + 29);
    const auto *sif_30 = buffer.data(sif + 30);
    const auto *sif_32 = buffer.data(sif + 32);
    const auto *sif_33 = buffer.data(sif + 33);
    const auto *sif_35 = buffer.data(sif + 35);
    const auto *sif_36 = buffer.data(sif + 36);
    const auto *sif_37 = buffer.data(sif + 37);
    const auto *sif_38 = buffer.data(sif + 38);
    const auto *sif_39 = buffer.data(sif + 39);
    const auto *sif_40 = buffer.data(sif + 40);
    const auto *sif_42 = buffer.data(sif + 42);
    const auto *sif_46 = buffer.data(sif + 46);
    const auto *sif_47 = buffer.data(sif + 47);
    const auto *sif_48 = buffer.data(sif + 48);
    const auto *sif_49 = buffer.data(sif + 49);
    const auto *sif_50 = buffer.data(sif + 50);
    const auto *sif_51 = buffer.data(sif + 51);
    const auto *sif_52 = buffer.data(sif + 52);
    const auto *sif_53 = buffer.data(sif + 53);
    const auto *sif_55 = buffer.data(sif + 55);
    const auto *sif_56 = buffer.data(sif + 56);
    const auto *sif_57 = buffer.data(sif + 57);
    const auto *sif_58 = buffer.data(sif + 58);
    const auto *sif_59 = buffer.data(sif + 59);
    const auto *sif_60 = buffer.data(sif + 60);
    const auto *sif_63 = buffer.data(sif + 63);
    const auto *sif_65 = buffer.data(sif + 65);
    const auto *sif_66 = buffer.data(sif + 66);
    const auto *sif_67 = buffer.data(sif + 67);
    const auto *sif_68 = buffer.data(sif + 68);
    const auto *sif_69 = buffer.data(sif + 69);
    const auto *sif_75 = buffer.data(sif + 75);
    const auto *sif_76 = buffer.data(sif + 76);
    const auto *sif_77 = buffer.data(sif + 77);
    const auto *sif_78 = buffer.data(sif + 78);
    const auto *sif_79 = buffer.data(sif + 79);

    const auto *sig1_0 = buffer.data(sig1 + 0);
    const auto *sig1_3 = buffer.data(sig1 + 3);
    const auto *sig1_5 = buffer.data(sig1 + 5);
    const auto *sig1_10 = buffer.data(sig1 + 10);
    const auto *sig1_14 = buffer.data(sig1 + 14);
    const auto *sig1_18 = buffer.data(sig1 + 18);
    const auto *sig1_25 = buffer.data(sig1 + 25);
    const auto *sig1_30 = buffer.data(sig1 + 30);
    const auto *sig1_35 = buffer.data(sig1 + 35);
    const auto *sig1_44 = buffer.data(sig1 + 44);
    const auto *sig1_45 = buffer.data(sig1 + 45);
    const auto *sig1_48 = buffer.data(sig1 + 48);
    const auto *sig1_55 = buffer.data(sig1 + 55);
    const auto *sig1_75 = buffer.data(sig1 + 75);
    const auto *sig1_78 = buffer.data(sig1 + 78);

    const auto *skd0_0 = buffer.data(skd0 + 0);
    const auto *skd0_3 = buffer.data(skd0 + 3);
    const auto *skd0_5 = buffer.data(skd0 + 5);
    const auto *skd0_9 = buffer.data(skd0 + 9);
    const auto *skd0_11 = buffer.data(skd0 + 11);
    const auto *skd0_17 = buffer.data(skd0 + 17);
    const auto *skd0_18 = buffer.data(skd0 + 18);
    const auto *skd0_21 = buffer.data(skd0 + 21);
    const auto *skd0_23 = buffer.data(skd0 + 23);
    const auto *skd0_29 = buffer.data(skd0 + 29);
    const auto *skd0_30 = buffer.data(skd0 + 30);
    const auto *skd0_33 = buffer.data(skd0 + 33);
    const auto *skd0_35 = buffer.data(skd0 + 35);
    const auto *skd0_36 = buffer.data(skd0 + 36);
    const auto *skd0_39 = buffer.data(skd0 + 39);
    const auto *skd0_41 = buffer.data(skd0 + 41);
    const auto *skd0_47 = buffer.data(skd0 + 47);

    const auto *skd1_0 = buffer.data(skd1 + 0);
    const auto *skd1_3 = buffer.data(skd1 + 3);
    const auto *skd1_5 = buffer.data(skd1 + 5);
    const auto *skd1_9 = buffer.data(skd1 + 9);
    const auto *skd1_11 = buffer.data(skd1 + 11);
    const auto *skd1_17 = buffer.data(skd1 + 17);
    const auto *skd1_18 = buffer.data(skd1 + 18);
    const auto *skd1_21 = buffer.data(skd1 + 21);
    const auto *skd1_23 = buffer.data(skd1 + 23);
    const auto *skd1_29 = buffer.data(skd1 + 29);
    const auto *skd1_30 = buffer.data(skd1 + 30);
    const auto *skd1_33 = buffer.data(skd1 + 33);
    const auto *skd1_35 = buffer.data(skd1 + 35);
    const auto *skd1_36 = buffer.data(skd1 + 36);
    const auto *skd1_39 = buffer.data(skd1 + 39);
    const auto *skd1_41 = buffer.data(skd1 + 41);
    const auto *skd1_47 = buffer.data(skd1 + 47);

    const auto *skf_0 = buffer.data(skf + 0);
    const auto *skf_2 = buffer.data(skf + 2);
    const auto *skf_3 = buffer.data(skf + 3);
    const auto *skf_5 = buffer.data(skf + 5);
    const auto *skf_6 = buffer.data(skf + 6);
    const auto *skf_7 = buffer.data(skf + 7);
    const auto *skf_8 = buffer.data(skf + 8);
    const auto *skf_9 = buffer.data(skf + 9);
    const auto *skf_10 = buffer.data(skf + 10);
    const auto *skf_12 = buffer.data(skf + 12);
    const auto *skf_16 = buffer.data(skf + 16);
    const auto *skf_17 = buffer.data(skf + 17);
    const auto *skf_18 = buffer.data(skf + 18);
    const auto *skf_19 = buffer.data(skf + 19);
    const auto *skf_20 = buffer.data(skf + 20);
    const auto *skf_22 = buffer.data(skf + 22);
    const auto *skf_26 = buffer.data(skf + 26);
    const auto *skf_27 = buffer.data(skf + 27);
    const auto *skf_28 = buffer.data(skf + 28);
    const auto *skf_29 = buffer.data(skf + 29);
    const auto *skf_30 = buffer.data(skf + 30);
    const auto *skf_32 = buffer.data(skf + 32);
    const auto *skf_33 = buffer.data(skf + 33);
    const auto *skf_35 = buffer.data(skf + 35);
    const auto *skf_36 = buffer.data(skf + 36);
    const auto *skf_37 = buffer.data(skf + 37);
    const auto *skf_38 = buffer.data(skf + 38);
    const auto *skf_39 = buffer.data(skf + 39);
    const auto *skf_40 = buffer.data(skf + 40);
    const auto *skf_42 = buffer.data(skf + 42);
    const auto *skf_46 = buffer.data(skf + 46);
    const auto *skf_47 = buffer.data(skf + 47);
    const auto *skf_48 = buffer.data(skf + 48);
    const auto *skf_49 = buffer.data(skf + 49);
    const auto *skf_50 = buffer.data(skf + 50);
    const auto *skf_52 = buffer.data(skf + 52);
    const auto *skf_53 = buffer.data(skf + 53);
    const auto *skf_55 = buffer.data(skf + 55);
    const auto *skf_56 = buffer.data(skf + 56);
    const auto *skf_57 = buffer.data(skf + 57);
    const auto *skf_58 = buffer.data(skf + 58);
    const auto *skf_59 = buffer.data(skf + 59);
    const auto *skf_60 = buffer.data(skf + 60);
    const auto *skf_62 = buffer.data(skf + 62);
    const auto *skf_63 = buffer.data(skf + 63);
    const auto *skf_65 = buffer.data(skf + 65);
    const auto *skf_66 = buffer.data(skf + 66);
    const auto *skf_67 = buffer.data(skf + 67);
    const auto *skf_68 = buffer.data(skf + 68);
    const auto *skf_69 = buffer.data(skf + 69);
    const auto *skf_70 = buffer.data(skf + 70);
    const auto *skf_72 = buffer.data(skf + 72);
    const auto *skf_75 = buffer.data(skf + 75);
    const auto *skf_76 = buffer.data(skf + 76);
    const auto *skf_77 = buffer.data(skf + 77);
    const auto *skf_78 = buffer.data(skf + 78);
    const auto *skf_79 = buffer.data(skf + 79);
    const auto *skf_80 = buffer.data(skf + 80);
    const auto *skf_82 = buffer.data(skf + 82);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, sif_0, sif_3, skd0_0, skd0_3, \
                         skd1_0, skd1_3, skf_0, skf_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sif_0[k]
                 + f_1 * skd0_0[k]
                 - f_2 * skd1_0[k]
                 + f_3 * pc_x[k] * skf_0[k];

        t_1[k] = f_3 * pc_y[k] * skf_0[k];

        t_2[k] = f_3 * pc_z[k] * skf_0[k];

        t_3[k] = f_0 * sif_3[k]
                 + f_4 * skd0_3[k]
                 - f_5 * skd1_3[k]
                 + f_3 * pc_x[k] * skf_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pc_x, pc_y, sif_5, sif_6, sif_7, skd0_5, skd1_5, \
                         skf_2, skf_5, skf_6, skf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * skf_2[k];

        t_5[k] = f_0 * sif_5[k]
                 + f_4 * skd0_5[k]
                 - f_5 * skd1_5[k]
                 + f_3 * pc_x[k] * skf_5[k];

        t_6[k] = f_0 * sif_6[k]
                 + f_3 * pc_x[k] * skf_6[k];

        t_7[k] = f_0 * sif_7[k]
                 + f_3 * pc_x[k] * skf_7[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pc_x, pc_y, pc_z, sif_8, sif_9, skd0_3, skd1_3, \
                         skf_6, skf_8, skf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_0 * sif_8[k]
                 + f_3 * pc_x[k] * skf_8[k];

        t_9[k] = f_0 * sif_9[k]
                 + f_3 * pc_x[k] * skf_9[k];

        t_10[k] = f_1 * skd0_3[k]
                  - f_2 * skd1_3[k]
                  + f_3 * pc_y[k] * skf_6[k];

        t_11[k] = f_3 * pc_z[k] * skf_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, pb_y, pc_y, pc_z, sig0_0, sif_0, \
                         sig1_0, skd0_5, skd1_5, skf_8, skf_9, skf_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_4 * skd0_5[k]
                  - f_5 * skd1_5[k]
                  + f_3 * pc_y[k] * skf_8[k];

        t_13[k] = f_3 * pc_y[k] * skf_9[k];

        t_14[k] = f_1 * skd0_5[k]
                  - f_2 * skd1_5[k]
                  + f_3 * pc_z[k] * skf_9[k];

        t_15[k] = pb_y[k] * sig0_0[k]
                  - f_6 * pc_y[k] * sig1_0[k];

        t_16[k] = f_7 * sif_0[k]
                  + f_3 * pc_y[k] * skf_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pb_y, pc_y, pc_z, sig0_3, sig0_5, sif_1, \
                         sif_2, sig1_3, sig1_5, skf_10, skf_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_3 * pc_z[k] * skf_10[k];

        t_18[k] = pb_y[k] * sig0_3[k]
                  + f_8 * sif_1[k]
                  - f_6 * pc_y[k] * sig1_3[k];

        t_19[k] = f_7 * sif_2[k]
                  + f_3 * pc_y[k] * skf_12[k];

        t_20[k] = pb_y[k] * sig0_5[k]
                  - f_6 * pc_y[k] * sig1_5[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_x, sif_16, sif_17, sif_18, sif_19, skf_16, \
                         skf_17, skf_18, skf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_9 * sif_16[k]
                  + f_3 * pc_x[k] * skf_16[k];

        t_22[k] = f_9 * sif_17[k]
                  + f_3 * pc_x[k] * skf_17[k];

        t_23[k] = f_9 * sif_18[k]
                  + f_3 * pc_x[k] * skf_18[k];

        t_24[k] = f_9 * sif_19[k]
                  + f_3 * pc_x[k] * skf_19[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pc_y, pc_z, sif_6, sif_8, sif_9, skd0_9, \
                         skd0_11, skd1_9, skd1_11, skf_16, skf_18, \
                         skf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_7 * sif_6[k]
                  + f_1 * skd0_9[k]
                  - f_2 * skd1_9[k]
                  + f_3 * pc_y[k] * skf_16[k];

        t_26[k] = f_3 * pc_z[k] * skf_16[k];

        t_27[k] = f_7 * sif_8[k]
                  + f_4 * skd0_11[k]
                  - f_5 * skd1_11[k]
                  + f_3 * pc_y[k] * skf_18[k];

        t_28[k] = f_7 * sif_9[k]
                  + f_3 * pc_y[k] * skf_19[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pb_y, pb_z, pc_y, pc_z, sig0_0, sig0_14, \
                         sif_0, sig1_0, sig1_14, skf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pb_y[k] * sig0_14[k]
                  - f_6 * pc_y[k] * sig1_14[k];

        t_30[k] = pb_z[k] * sig0_0[k]
                  - f_6 * pc_z[k] * sig1_0[k];

        t_31[k] = f_3 * pc_y[k] * skf_20[k];

        t_32[k] = f_7 * sif_0[k]
                  + f_3 * pc_z[k] * skf_20[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pb_z, pc_x, pc_y, pc_z, sig0_3, sig0_5, \
                         sif_2, sif_26, sig1_3, sig1_5, skf_22, \
                         skf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pb_z[k] * sig0_3[k]
                  - f_6 * pc_z[k] * sig1_3[k];

        t_34[k] = f_3 * pc_y[k] * skf_22[k];

        t_35[k] = pb_z[k] * sig0_5[k]
                  + f_8 * sif_2[k]
                  - f_6 * pc_z[k] * sig1_5[k];

        t_36[k] = f_9 * sif_26[k]
                  + f_3 * pc_x[k] * skf_26[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pb_z, pc_x, pc_z, sig0_10, sif_27, sif_28, \
                         sif_29, sig1_10, skf_27, skf_28, skf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_9 * sif_27[k]
                  + f_3 * pc_x[k] * skf_27[k];

        t_38[k] = f_9 * sif_28[k]
                  + f_3 * pc_x[k] * skf_28[k];

        t_39[k] = f_9 * sif_29[k]
                  + f_3 * pc_x[k] * skf_29[k];

        t_40[k] = pb_z[k] * sig0_10[k]
                  - f_6 * pc_z[k] * sig1_10[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, pc_y, pc_z, sif_6, sif_9, skd0_17, skd1_17, \
                         skf_26, skf_28, skf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_7 * sif_6[k]
                  + f_3 * pc_z[k] * skf_26[k];

        t_42[k] = f_4 * skd0_17[k]
                  - f_5 * skd1_17[k]
                  + f_3 * pc_y[k] * skf_28[k];

        t_43[k] = f_3 * pc_y[k] * skf_29[k];

        t_44[k] = f_7 * sif_9[k]
                  + f_1 * skd0_17[k]
                  - f_2 * skd1_17[k]
                  + f_3 * pc_z[k] * skf_29[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pc_x, pc_y, pc_z, sif_10, sif_30, sif_33, \
                         skd0_18, skd0_21, skd1_18, skd1_21, skf_30, \
                         skf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_10 * sif_30[k]
                  + f_1 * skd0_18[k]
                  - f_2 * skd1_18[k]
                  + f_3 * pc_x[k] * skf_30[k];

        t_46[k] = f_8 * sif_10[k]
                  + f_3 * pc_y[k] * skf_30[k];

        t_47[k] = f_3 * pc_z[k] * skf_30[k];

        t_48[k] = f_10 * sif_33[k]
                  + f_4 * skd0_21[k]
                  - f_5 * skd1_21[k]
                  + f_3 * pc_x[k] * skf_33[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pc_x, pc_y, sif_12, sif_35, sif_36, sif_37, \
                         skd0_23, skd1_23, skf_32, skf_35, skf_36, \
                         skf_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_8 * sif_12[k]
                  + f_3 * pc_y[k] * skf_32[k];

        t_50[k] = f_10 * sif_35[k]
                  + f_4 * skd0_23[k]
                  - f_5 * skd1_23[k]
                  + f_3 * pc_x[k] * skf_35[k];

        t_51[k] = f_10 * sif_36[k]
                  + f_3 * pc_x[k] * skf_36[k];

        t_52[k] = f_10 * sif_37[k]
                  + f_3 * pc_x[k] * skf_37[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pc_x, pc_y, pc_z, sif_16, sif_38, sif_39, \
                         skd0_21, skd1_21, skf_36, skf_38, skf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_10 * sif_38[k]
                  + f_3 * pc_x[k] * skf_38[k];

        t_54[k] = f_10 * sif_39[k]
                  + f_3 * pc_x[k] * skf_39[k];

        t_55[k] = f_8 * sif_16[k]
                  + f_1 * skd0_21[k]
                  - f_2 * skd1_21[k]
                  + f_3 * pc_y[k] * skf_36[k];

        t_56[k] = f_3 * pc_z[k] * skf_36[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pb_y, pc_y, pc_z, sig0_30, sif_18, sif_19, \
                         sig1_30, skd0_23, skd1_23, skf_38, skf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_8 * sif_18[k]
                  + f_4 * skd0_23[k]
                  - f_5 * skd1_23[k]
                  + f_3 * pc_y[k] * skf_38[k];

        t_58[k] = f_8 * sif_19[k]
                  + f_3 * pc_y[k] * skf_39[k];

        t_59[k] = f_1 * skd0_23[k]
                  - f_2 * skd1_23[k]
                  + f_3 * pc_z[k] * skf_39[k];

        t_60[k] = pb_y[k] * sig0_30[k]
                  - f_6 * pc_y[k] * sig1_30[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pb_z, pc_y, pc_z, sig0_18, sif_10, sif_20, \
                         sif_22, sig1_18, skf_40, skf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_7 * sif_20[k]
                  + f_3 * pc_y[k] * skf_40[k];

        t_62[k] = f_7 * sif_10[k]
                  + f_3 * pc_z[k] * skf_40[k];

        t_63[k] = pb_z[k] * sig0_18[k]
                  - f_6 * pc_z[k] * sig1_18[k];

        t_64[k] = f_7 * sif_22[k]
                  + f_3 * pc_y[k] * skf_42[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pb_y, pc_x, pc_y, sig0_35, sif_46, sif_47, \
                         sif_48, sig1_35, skf_46, skf_47, skf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pb_y[k] * sig0_35[k]
                  - f_6 * pc_y[k] * sig1_35[k];

        t_66[k] = f_10 * sif_46[k]
                  + f_3 * pc_x[k] * skf_46[k];

        t_67[k] = f_10 * sif_47[k]
                  + f_3 * pc_x[k] * skf_47[k];

        t_68[k] = f_10 * sif_48[k]
                  + f_3 * pc_x[k] * skf_48[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pb_z, pc_x, pc_z, sig0_25, sif_16, sif_49, sig1_25, \
                         skf_46, skf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_10 * sif_49[k]
                  + f_3 * pc_x[k] * skf_49[k];

        t_70[k] = pb_z[k] * sig0_25[k]
                  - f_6 * pc_z[k] * sig1_25[k];

        t_71[k] = f_7 * sif_16[k]
                  + f_3 * pc_z[k] * skf_46[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pb_y, pc_y, sig0_44, sif_28, sif_29, sig1_44, \
                         skd0_29, skd1_29, skf_48, skf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_7 * sif_28[k]
                  + f_4 * skd0_29[k]
                  - f_5 * skd1_29[k]
                  + f_3 * pc_y[k] * skf_48[k];

        t_73[k] = f_7 * sif_29[k]
                  + f_3 * pc_y[k] * skf_49[k];

        t_74[k] = pb_y[k] * sig0_44[k]
                  - f_6 * pc_y[k] * sig1_44[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, pc_x, pc_y, pc_z, sif_20, sif_50, sif_53, \
                         skd0_30, skd0_33, skd1_30, skd1_33, skf_50, \
                         skf_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_10 * sif_50[k]
                  + f_1 * skd0_30[k]
                  - f_2 * skd1_30[k]
                  + f_3 * pc_x[k] * skf_50[k];

        t_76[k] = f_3 * pc_y[k] * skf_50[k];

        t_77[k] = f_8 * sif_20[k]
                  + f_3 * pc_z[k] * skf_50[k];

        t_78[k] = f_10 * sif_53[k]
                  + f_4 * skd0_33[k]
                  - f_5 * skd1_33[k]
                  + f_3 * pc_x[k] * skf_53[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, t_82, pc_x, pc_y, sif_55, sif_56, sif_57, skd0_35, \
                         skd1_35, skf_52, skf_55, skf_56, skf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_3 * pc_y[k] * skf_52[k];

        t_80[k] = f_10 * sif_55[k]
                  + f_4 * skd0_35[k]
                  - f_5 * skd1_35[k]
                  + f_3 * pc_x[k] * skf_55[k];

        t_81[k] = f_10 * sif_56[k]
                  + f_3 * pc_x[k] * skf_56[k];

        t_82[k] = f_10 * sif_57[k]
                  + f_3 * pc_x[k] * skf_57[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pc_x, pc_y, pc_z, sif_26, sif_58, sif_59, \
                         skd0_33, skd1_33, skf_56, skf_58, skf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_10 * sif_58[k]
                  + f_3 * pc_x[k] * skf_58[k];

        t_84[k] = f_10 * sif_59[k]
                  + f_3 * pc_x[k] * skf_59[k];

        t_85[k] = f_1 * skd0_33[k]
                  - f_2 * skd1_33[k]
                  + f_3 * pc_y[k] * skf_56[k];

        t_86[k] = f_8 * sif_26[k]
                  + f_3 * pc_z[k] * skf_56[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pc_x, pc_y, pc_z, sif_29, sif_60, skd0_35, \
                         skd0_36, skd1_35, skd1_36, skf_58, skf_59, \
                         skf_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_4 * skd0_35[k]
                  - f_5 * skd1_35[k]
                  + f_3 * pc_y[k] * skf_58[k];

        t_88[k] = f_3 * pc_y[k] * skf_59[k];

        t_89[k] = f_8 * sif_29[k]
                  + f_1 * skd0_35[k]
                  - f_2 * skd1_35[k]
                  + f_3 * pc_z[k] * skf_59[k];

        t_90[k] = f_11 * sif_60[k]
                  + f_1 * skd0_36[k]
                  - f_2 * skd1_36[k]
                  + f_3 * pc_x[k] * skf_60[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pc_x, pc_y, pc_z, sif_30, sif_32, sif_63, \
                         skd0_39, skd1_39, skf_60, skf_62, skf_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_12 * sif_30[k]
                  + f_3 * pc_y[k] * skf_60[k];

        t_92[k] = f_3 * pc_z[k] * skf_60[k];

        t_93[k] = f_11 * sif_63[k]
                  + f_4 * skd0_39[k]
                  - f_5 * skd1_39[k]
                  + f_3 * pc_x[k] * skf_63[k];

        t_94[k] = f_12 * sif_32[k]
                  + f_3 * pc_y[k] * skf_62[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pc_x, sif_65, sif_66, sif_67, sif_68, \
                         skd0_41, skd1_41, skf_65, skf_66, skf_67, \
                         skf_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_11 * sif_65[k]
                  + f_4 * skd0_41[k]
                  - f_5 * skd1_41[k]
                  + f_3 * pc_x[k] * skf_65[k];

        t_96[k] = f_11 * sif_66[k]
                  + f_3 * pc_x[k] * skf_66[k];

        t_97[k] = f_11 * sif_67[k]
                  + f_3 * pc_x[k] * skf_67[k];

        t_98[k] = f_11 * sif_68[k]
                  + f_3 * pc_x[k] * skf_68[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, pc_x, pc_y, pc_z, sif_36, sif_69, skd0_39, \
                         skd1_39, skf_66, skf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_11 * sif_69[k]
                  + f_3 * pc_x[k] * skf_69[k];

        t_100[k] = f_12 * sif_36[k]
                   + f_1 * skd0_39[k]
                   - f_2 * skd1_39[k]
                   + f_3 * pc_y[k] * skf_66[k];

        t_101[k] = f_3 * pc_z[k] * skf_66[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, pb_z, pc_y, pc_z, sig0_45, sif_38, \
                         sif_39, sig1_45, skd0_41, skd1_41, skf_68, \
                         skf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = f_12 * sif_38[k]
                   + f_4 * skd0_41[k]
                   - f_5 * skd1_41[k]
                   + f_3 * pc_y[k] * skf_68[k];

        t_103[k] = f_12 * sif_39[k]
                   + f_3 * pc_y[k] * skf_69[k];

        t_104[k] = f_1 * skd0_41[k]
                   - f_2 * skd1_41[k]
                   + f_3 * pc_z[k] * skf_69[k];

        t_105[k] = pb_z[k] * sig0_45[k]
                   - f_6 * pc_z[k] * sig1_45[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pb_z, pc_y, pc_z, sig0_48, sif_30, \
                         sif_40, sif_42, sig1_48, skf_70, skf_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_8 * sif_40[k]
                   + f_3 * pc_y[k] * skf_70[k];

        t_107[k] = f_7 * sif_30[k]
                   + f_3 * pc_z[k] * skf_70[k];

        t_108[k] = pb_z[k] * sig0_48[k]
                   - f_6 * pc_z[k] * sig1_48[k];

        t_109[k] = f_8 * sif_42[k]
                   + f_3 * pc_y[k] * skf_72[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, pc_x, sif_75, sif_76, sif_77, sif_78, \
                         skd0_47, skd1_47, skf_75, skf_76, skf_77, \
                         skf_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_11 * sif_75[k]
                   + f_4 * skd0_47[k]
                   - f_5 * skd1_47[k]
                   + f_3 * pc_x[k] * skf_75[k];

        t_111[k] = f_11 * sif_76[k]
                   + f_3 * pc_x[k] * skf_76[k];

        t_112[k] = f_11 * sif_77[k]
                   + f_3 * pc_x[k] * skf_77[k];

        t_113[k] = f_11 * sif_78[k]
                   + f_3 * pc_x[k] * skf_78[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, pb_z, pc_x, pc_z, sig0_55, sif_36, sif_79, \
                         sig1_55, skf_76, skf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_11 * sif_79[k]
                   + f_3 * pc_x[k] * skf_79[k];

        t_115[k] = pb_z[k] * sig0_55[k]
                   - f_6 * pc_z[k] * sig1_55[k];

        t_116[k] = f_7 * sif_36[k]
                   + f_3 * pc_z[k] * skf_76[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pb_y, pc_y, pc_z, sig0_75, sif_39, \
                         sif_48, sif_49, sig1_75, skd0_47, skd1_47, skf_78, \
                         skf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_8 * sif_48[k]
                   + f_4 * skd0_47[k]
                   - f_5 * skd1_47[k]
                   + f_3 * pc_y[k] * skf_78[k];

        t_118[k] = f_8 * sif_49[k]
                   + f_3 * pc_y[k] * skf_79[k];

        t_119[k] = f_7 * sif_39[k]
                   + f_1 * skd0_47[k]
                   - f_2 * skd1_47[k]
                   + f_3 * pc_z[k] * skf_79[k];

        t_120[k] = pb_y[k] * sig0_75[k]
                   - f_6 * pc_y[k] * sig1_75[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pb_y, pc_y, pc_z, sig0_78, sif_40, \
                         sif_50, sif_51, sif_52, sig1_78, skf_80, \
                         skf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_7 * sif_50[k]
                   + f_3 * pc_y[k] * skf_80[k];

        t_122[k] = f_8 * sif_40[k]
                   + f_3 * pc_z[k] * skf_80[k];

        t_123[k] = pb_y[k] * sig0_78[k]
                   + f_8 * sif_51[k]
                   - f_6 * pc_y[k] * sig1_78[k];

        t_124[k] = f_7 * sif_52[k]
                   + f_3 * pc_y[k] * skf_82[k];
    }
}

static auto
compute_prim_skg_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sig0,
                                                          const size_t sif, const size_t sig1,
                                                          const size_t skd0, const size_t skd1,
                                                          const size_t skf, const size_t ncols,
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
    const auto f_10 = 2.5 / q;
    const auto f_11 = 2.0 / q;
    const auto f_12 = 1.5 / q;

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

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sig0_80 = buffer.data(sig0 + 80);
    const auto *sig0_89 = buffer.data(sig0 + 89);
    const auto *sig0_90 = buffer.data(sig0 + 90);
    const auto *sig0_93 = buffer.data(sig0 + 93);
    const auto *sig0_100 = buffer.data(sig0 + 100);
    const auto *sig0_135 = buffer.data(sig0 + 135);
    const auto *sig0_138 = buffer.data(sig0 + 138);
    const auto *sig0_140 = buffer.data(sig0 + 140);
    const auto *sig0_149 = buffer.data(sig0 + 149);
    const auto *sig0_150 = buffer.data(sig0 + 150);
    const auto *sig0_153 = buffer.data(sig0 + 153);

    const auto *sif_46 = buffer.data(sif + 46);
    const auto *sif_50 = buffer.data(sif + 50);
    const auto *sif_56 = buffer.data(sif + 56);
    const auto *sif_58 = buffer.data(sif + 58);
    const auto *sif_59 = buffer.data(sif + 59);
    const auto *sif_60 = buffer.data(sif + 60);
    const auto *sif_62 = buffer.data(sif + 62);
    const auto *sif_66 = buffer.data(sif + 66);
    const auto *sif_68 = buffer.data(sif + 68);
    const auto *sif_69 = buffer.data(sif + 69);
    const auto *sif_70 = buffer.data(sif + 70);
    const auto *sif_72 = buffer.data(sif + 72);
    const auto *sif_76 = buffer.data(sif + 76);
    const auto *sif_78 = buffer.data(sif + 78);
    const auto *sif_79 = buffer.data(sif + 79);
    const auto *sif_80 = buffer.data(sif + 80);
    const auto *sif_82 = buffer.data(sif + 82);
    const auto *sif_86 = buffer.data(sif + 86);
    const auto *sif_87 = buffer.data(sif + 87);
    const auto *sif_88 = buffer.data(sif + 88);
    const auto *sif_89 = buffer.data(sif + 89);
    const auto *sif_90 = buffer.data(sif + 90);
    const auto *sif_91 = buffer.data(sif + 91);
    const auto *sif_92 = buffer.data(sif + 92);
    const auto *sif_93 = buffer.data(sif + 93);
    const auto *sif_95 = buffer.data(sif + 95);
    const auto *sif_96 = buffer.data(sif + 96);
    const auto *sif_97 = buffer.data(sif + 97);
    const auto *sif_98 = buffer.data(sif + 98);
    const auto *sif_99 = buffer.data(sif + 99);
    const auto *sif_100 = buffer.data(sif + 100);
    const auto *sif_102 = buffer.data(sif + 102);
    const auto *sif_103 = buffer.data(sif + 103);
    const auto *sif_105 = buffer.data(sif + 105);
    const auto *sif_106 = buffer.data(sif + 106);
    const auto *sif_107 = buffer.data(sif + 107);
    const auto *sif_108 = buffer.data(sif + 108);
    const auto *sif_109 = buffer.data(sif + 109);
    const auto *sif_110 = buffer.data(sif + 110);
    const auto *sif_112 = buffer.data(sif + 112);
    const auto *sif_115 = buffer.data(sif + 115);
    const auto *sif_116 = buffer.data(sif + 116);
    const auto *sif_117 = buffer.data(sif + 117);
    const auto *sif_118 = buffer.data(sif + 118);
    const auto *sif_119 = buffer.data(sif + 119);
    const auto *sif_120 = buffer.data(sif + 120);
    const auto *sif_123 = buffer.data(sif + 123);
    const auto *sif_125 = buffer.data(sif + 125);
    const auto *sif_126 = buffer.data(sif + 126);
    const auto *sif_127 = buffer.data(sif + 127);
    const auto *sif_128 = buffer.data(sif + 128);
    const auto *sif_129 = buffer.data(sif + 129);
    const auto *sif_136 = buffer.data(sif + 136);
    const auto *sif_137 = buffer.data(sif + 137);
    const auto *sif_138 = buffer.data(sif + 138);
    const auto *sif_139 = buffer.data(sif + 139);
    const auto *sif_140 = buffer.data(sif + 140);
    const auto *sif_143 = buffer.data(sif + 143);
    const auto *sif_145 = buffer.data(sif + 145);
    const auto *sif_146 = buffer.data(sif + 146);
    const auto *sif_147 = buffer.data(sif + 147);
    const auto *sif_148 = buffer.data(sif + 148);
    const auto *sif_149 = buffer.data(sif + 149);
    const auto *sif_150 = buffer.data(sif + 150);
    const auto *sif_153 = buffer.data(sif + 153);
    const auto *sif_155 = buffer.data(sif + 155);
    const auto *sif_156 = buffer.data(sif + 156);
    const auto *sif_157 = buffer.data(sif + 157);
    const auto *sif_158 = buffer.data(sif + 158);
    const auto *sif_159 = buffer.data(sif + 159);

    const auto *sig1_80 = buffer.data(sig1 + 80);
    const auto *sig1_89 = buffer.data(sig1 + 89);
    const auto *sig1_90 = buffer.data(sig1 + 90);
    const auto *sig1_93 = buffer.data(sig1 + 93);
    const auto *sig1_100 = buffer.data(sig1 + 100);
    const auto *sig1_135 = buffer.data(sig1 + 135);
    const auto *sig1_138 = buffer.data(sig1 + 138);
    const auto *sig1_140 = buffer.data(sig1 + 140);
    const auto *sig1_149 = buffer.data(sig1 + 149);
    const auto *sig1_150 = buffer.data(sig1 + 150);
    const auto *sig1_153 = buffer.data(sig1 + 153);

    const auto *skd0_51 = buffer.data(skd0 + 51);
    const auto *skd0_53 = buffer.data(skd0 + 53);
    const auto *skd0_54 = buffer.data(skd0 + 54);
    const auto *skd0_57 = buffer.data(skd0 + 57);
    const auto *skd0_59 = buffer.data(skd0 + 59);
    const auto *skd0_60 = buffer.data(skd0 + 60);
    const auto *skd0_63 = buffer.data(skd0 + 63);
    const auto *skd0_65 = buffer.data(skd0 + 65);
    const auto *skd0_71 = buffer.data(skd0 + 71);
    const auto *skd0_72 = buffer.data(skd0 + 72);
    const auto *skd0_75 = buffer.data(skd0 + 75);
    const auto *skd0_77 = buffer.data(skd0 + 77);
    const auto *skd0_81 = buffer.data(skd0 + 81);
    const auto *skd0_83 = buffer.data(skd0 + 83);
    const auto *skd0_84 = buffer.data(skd0 + 84);
    const auto *skd0_87 = buffer.data(skd0 + 87);
    const auto *skd0_89 = buffer.data(skd0 + 89);
    const auto *skd0_90 = buffer.data(skd0 + 90);
    const auto *skd0_93 = buffer.data(skd0 + 93);
    const auto *skd0_95 = buffer.data(skd0 + 95);

    const auto *skd1_51 = buffer.data(skd1 + 51);
    const auto *skd1_53 = buffer.data(skd1 + 53);
    const auto *skd1_54 = buffer.data(skd1 + 54);
    const auto *skd1_57 = buffer.data(skd1 + 57);
    const auto *skd1_59 = buffer.data(skd1 + 59);
    const auto *skd1_60 = buffer.data(skd1 + 60);
    const auto *skd1_63 = buffer.data(skd1 + 63);
    const auto *skd1_65 = buffer.data(skd1 + 65);
    const auto *skd1_71 = buffer.data(skd1 + 71);
    const auto *skd1_72 = buffer.data(skd1 + 72);
    const auto *skd1_75 = buffer.data(skd1 + 75);
    const auto *skd1_77 = buffer.data(skd1 + 77);
    const auto *skd1_81 = buffer.data(skd1 + 81);
    const auto *skd1_83 = buffer.data(skd1 + 83);
    const auto *skd1_84 = buffer.data(skd1 + 84);
    const auto *skd1_87 = buffer.data(skd1 + 87);
    const auto *skd1_89 = buffer.data(skd1 + 89);
    const auto *skd1_90 = buffer.data(skd1 + 90);
    const auto *skd1_93 = buffer.data(skd1 + 93);
    const auto *skd1_95 = buffer.data(skd1 + 95);

    const auto *skf_86 = buffer.data(skf + 86);
    const auto *skf_87 = buffer.data(skf + 87);
    const auto *skf_88 = buffer.data(skf + 88);
    const auto *skf_89 = buffer.data(skf + 89);
    const auto *skf_90 = buffer.data(skf + 90);
    const auto *skf_92 = buffer.data(skf + 92);
    const auto *skf_93 = buffer.data(skf + 93);
    const auto *skf_95 = buffer.data(skf + 95);
    const auto *skf_96 = buffer.data(skf + 96);
    const auto *skf_97 = buffer.data(skf + 97);
    const auto *skf_98 = buffer.data(skf + 98);
    const auto *skf_99 = buffer.data(skf + 99);
    const auto *skf_100 = buffer.data(skf + 100);
    const auto *skf_102 = buffer.data(skf + 102);
    const auto *skf_103 = buffer.data(skf + 103);
    const auto *skf_105 = buffer.data(skf + 105);
    const auto *skf_106 = buffer.data(skf + 106);
    const auto *skf_107 = buffer.data(skf + 107);
    const auto *skf_108 = buffer.data(skf + 108);
    const auto *skf_109 = buffer.data(skf + 109);
    const auto *skf_110 = buffer.data(skf + 110);
    const auto *skf_112 = buffer.data(skf + 112);
    const auto *skf_115 = buffer.data(skf + 115);
    const auto *skf_116 = buffer.data(skf + 116);
    const auto *skf_117 = buffer.data(skf + 117);
    const auto *skf_118 = buffer.data(skf + 118);
    const auto *skf_119 = buffer.data(skf + 119);
    const auto *skf_120 = buffer.data(skf + 120);
    const auto *skf_122 = buffer.data(skf + 122);
    const auto *skf_123 = buffer.data(skf + 123);
    const auto *skf_125 = buffer.data(skf + 125);
    const auto *skf_126 = buffer.data(skf + 126);
    const auto *skf_127 = buffer.data(skf + 127);
    const auto *skf_128 = buffer.data(skf + 128);
    const auto *skf_129 = buffer.data(skf + 129);
    const auto *skf_130 = buffer.data(skf + 130);
    const auto *skf_132 = buffer.data(skf + 132);
    const auto *skf_136 = buffer.data(skf + 136);
    const auto *skf_137 = buffer.data(skf + 137);
    const auto *skf_138 = buffer.data(skf + 138);
    const auto *skf_139 = buffer.data(skf + 139);
    const auto *skf_140 = buffer.data(skf + 140);
    const auto *skf_142 = buffer.data(skf + 142);
    const auto *skf_143 = buffer.data(skf + 143);
    const auto *skf_145 = buffer.data(skf + 145);
    const auto *skf_146 = buffer.data(skf + 146);
    const auto *skf_147 = buffer.data(skf + 147);
    const auto *skf_148 = buffer.data(skf + 148);
    const auto *skf_149 = buffer.data(skf + 149);
    const auto *skf_150 = buffer.data(skf + 150);
    const auto *skf_152 = buffer.data(skf + 152);
    const auto *skf_153 = buffer.data(skf + 153);
    const auto *skf_155 = buffer.data(skf + 155);
    const auto *skf_156 = buffer.data(skf + 156);
    const auto *skf_157 = buffer.data(skf + 157);
    const auto *skf_158 = buffer.data(skf + 158);
    const auto *skf_159 = buffer.data(skf + 159);
    const auto *skf_160 = buffer.data(skf + 160);
    const auto *skf_162 = buffer.data(skf + 162);

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pb_y, pc_x, pc_y, sig0_80, sif_86, \
                         sif_87, sif_88, sig1_80, skf_86, skf_87, \
                         skf_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = pb_y[k] * sig0_80[k]
                   - f_6 * pc_y[k] * sig1_80[k];

        t_126[k] = f_11 * sif_86[k]
                   + f_3 * pc_x[k] * skf_86[k];

        t_127[k] = f_11 * sif_87[k]
                   + f_3 * pc_x[k] * skf_87[k];

        t_128[k] = f_11 * sif_88[k]
                   + f_3 * pc_x[k] * skf_88[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, pc_x, pc_y, pc_z, sif_46, sif_56, sif_89, \
                         skd0_51, skd1_51, skf_86, skf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_11 * sif_89[k]
                   + f_3 * pc_x[k] * skf_89[k];

        t_130[k] = f_7 * sif_56[k]
                   + f_1 * skd0_51[k]
                   - f_2 * skd1_51[k]
                   + f_3 * pc_y[k] * skf_86[k];

        t_131[k] = f_8 * sif_46[k]
                   + f_3 * pc_z[k] * skf_86[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, pb_y, pc_y, sig0_89, sif_58, sif_59, sig1_89, \
                         skd0_53, skd1_53, skf_88, skf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_7 * sif_58[k]
                   + f_4 * skd0_53[k]
                   - f_5 * skd1_53[k]
                   + f_3 * pc_y[k] * skf_88[k];

        t_133[k] = f_7 * sif_59[k]
                   + f_3 * pc_y[k] * skf_89[k];

        t_134[k] = pb_y[k] * sig0_89[k]
                   - f_6 * pc_y[k] * sig1_89[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, pc_x, pc_y, pc_z, sif_50, sif_90, sif_93, \
                         skd0_54, skd0_57, skd1_54, skd1_57, skf_90, \
                         skf_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_11 * sif_90[k]
                   + f_1 * skd0_54[k]
                   - f_2 * skd1_54[k]
                   + f_3 * pc_x[k] * skf_90[k];

        t_136[k] = f_3 * pc_y[k] * skf_90[k];

        t_137[k] = f_12 * sif_50[k]
                   + f_3 * pc_z[k] * skf_90[k];

        t_138[k] = f_11 * sif_93[k]
                   + f_4 * skd0_57[k]
                   - f_5 * skd1_57[k]
                   + f_3 * pc_x[k] * skf_93[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pc_x, pc_y, sif_95, sif_96, sif_97, \
                         skd0_59, skd1_59, skf_92, skf_95, skf_96, \
                         skf_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_3 * pc_y[k] * skf_92[k];

        t_140[k] = f_11 * sif_95[k]
                   + f_4 * skd0_59[k]
                   - f_5 * skd1_59[k]
                   + f_3 * pc_x[k] * skf_95[k];

        t_141[k] = f_11 * sif_96[k]
                   + f_3 * pc_x[k] * skf_96[k];

        t_142[k] = f_11 * sif_97[k]
                   + f_3 * pc_x[k] * skf_97[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, t_146, pc_x, pc_y, pc_z, sif_56, sif_98, sif_99, \
                         skd0_57, skd1_57, skf_96, skf_98, skf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_11 * sif_98[k]
                   + f_3 * pc_x[k] * skf_98[k];

        t_144[k] = f_11 * sif_99[k]
                   + f_3 * pc_x[k] * skf_99[k];

        t_145[k] = f_1 * skd0_57[k]
                   - f_2 * skd1_57[k]
                   + f_3 * pc_y[k] * skf_96[k];

        t_146[k] = f_12 * sif_56[k]
                   + f_3 * pc_z[k] * skf_96[k];
    }

#pragma omp simd aligned(t_147, t_148, t_149, t_150, pc_x, pc_y, pc_z, sif_59, sif_100, \
                         skd0_59, skd0_60, skd1_59, skd1_60, skf_98, skf_99, \
                         skf_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_147[k] = f_4 * skd0_59[k]
                   - f_5 * skd1_59[k]
                   + f_3 * pc_y[k] * skf_98[k];

        t_148[k] = f_3 * pc_y[k] * skf_99[k];

        t_149[k] = f_12 * sif_59[k]
                   + f_1 * skd0_59[k]
                   - f_2 * skd1_59[k]
                   + f_3 * pc_z[k] * skf_99[k];

        t_150[k] = f_12 * sif_100[k]
                   + f_1 * skd0_60[k]
                   - f_2 * skd1_60[k]
                   + f_3 * pc_x[k] * skf_100[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pc_x, pc_y, pc_z, sif_60, sif_62, \
                         sif_103, skd0_63, skd1_63, skf_100, skf_102, \
                         skf_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_11 * sif_60[k]
                   + f_3 * pc_y[k] * skf_100[k];

        t_152[k] = f_3 * pc_z[k] * skf_100[k];

        t_153[k] = f_12 * sif_103[k]
                   + f_4 * skd0_63[k]
                   - f_5 * skd1_63[k]
                   + f_3 * pc_x[k] * skf_103[k];

        t_154[k] = f_11 * sif_62[k]
                   + f_3 * pc_y[k] * skf_102[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pc_x, sif_105, sif_106, sif_107, sif_108, \
                         skd0_65, skd1_65, skf_105, skf_106, skf_107, \
                         skf_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_12 * sif_105[k]
                   + f_4 * skd0_65[k]
                   - f_5 * skd1_65[k]
                   + f_3 * pc_x[k] * skf_105[k];

        t_156[k] = f_12 * sif_106[k]
                   + f_3 * pc_x[k] * skf_106[k];

        t_157[k] = f_12 * sif_107[k]
                   + f_3 * pc_x[k] * skf_107[k];

        t_158[k] = f_12 * sif_108[k]
                   + f_3 * pc_x[k] * skf_108[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, pc_x, pc_y, pc_z, sif_66, sif_109, skd0_63, \
                         skd1_63, skf_106, skf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_12 * sif_109[k]
                   + f_3 * pc_x[k] * skf_109[k];

        t_160[k] = f_11 * sif_66[k]
                   + f_1 * skd0_63[k]
                   - f_2 * skd1_63[k]
                   + f_3 * pc_y[k] * skf_106[k];

        t_161[k] = f_3 * pc_z[k] * skf_106[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pb_z, pc_y, pc_z, sig0_90, sif_68, \
                         sif_69, sig1_90, skd0_65, skd1_65, skf_108, \
                         skf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_11 * sif_68[k]
                   + f_4 * skd0_65[k]
                   - f_5 * skd1_65[k]
                   + f_3 * pc_y[k] * skf_108[k];

        t_163[k] = f_11 * sif_69[k]
                   + f_3 * pc_y[k] * skf_109[k];

        t_164[k] = f_1 * skd0_65[k]
                   - f_2 * skd1_65[k]
                   + f_3 * pc_z[k] * skf_109[k];

        t_165[k] = pb_z[k] * sig0_90[k]
                   - f_6 * pc_z[k] * sig1_90[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pb_z, pc_y, pc_z, sig0_93, sif_60, \
                         sif_70, sif_72, sig1_93, skf_110, skf_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_12 * sif_70[k]
                   + f_3 * pc_y[k] * skf_110[k];

        t_167[k] = f_7 * sif_60[k]
                   + f_3 * pc_z[k] * skf_110[k];

        t_168[k] = pb_z[k] * sig0_93[k]
                   - f_6 * pc_z[k] * sig1_93[k];

        t_169[k] = f_12 * sif_72[k]
                   + f_3 * pc_y[k] * skf_112[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pc_x, sif_115, sif_116, sif_117, sif_118, \
                         skd0_71, skd1_71, skf_115, skf_116, skf_117, \
                         skf_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_12 * sif_115[k]
                   + f_4 * skd0_71[k]
                   - f_5 * skd1_71[k]
                   + f_3 * pc_x[k] * skf_115[k];

        t_171[k] = f_12 * sif_116[k]
                   + f_3 * pc_x[k] * skf_116[k];

        t_172[k] = f_12 * sif_117[k]
                   + f_3 * pc_x[k] * skf_117[k];

        t_173[k] = f_12 * sif_118[k]
                   + f_3 * pc_x[k] * skf_118[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pb_z, pc_x, pc_z, sig0_100, sif_66, sif_119, \
                         sig1_100, skf_116, skf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_12 * sif_119[k]
                   + f_3 * pc_x[k] * skf_119[k];

        t_175[k] = pb_z[k] * sig0_100[k]
                   - f_6 * pc_z[k] * sig1_100[k];

        t_176[k] = f_7 * sif_66[k]
                   + f_3 * pc_z[k] * skf_116[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pc_y, pc_z, sif_69, sif_78, sif_79, skd0_71, \
                         skd1_71, skf_118, skf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_12 * sif_78[k]
                   + f_4 * skd0_71[k]
                   - f_5 * skd1_71[k]
                   + f_3 * pc_y[k] * skf_118[k];

        t_178[k] = f_12 * sif_79[k]
                   + f_3 * pc_y[k] * skf_119[k];

        t_179[k] = f_7 * sif_69[k]
                   + f_1 * skd0_71[k]
                   - f_2 * skd1_71[k]
                   + f_3 * pc_z[k] * skf_119[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, pc_x, pc_y, pc_z, sif_70, sif_80, sif_120, \
                         skd0_72, skd1_72, skf_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_12 * sif_120[k]
                   + f_1 * skd0_72[k]
                   - f_2 * skd1_72[k]
                   + f_3 * pc_x[k] * skf_120[k];

        t_181[k] = f_8 * sif_80[k]
                   + f_3 * pc_y[k] * skf_120[k];

        t_182[k] = f_8 * sif_70[k]
                   + f_3 * pc_z[k] * skf_120[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pc_x, pc_y, sif_82, sif_123, sif_125, skd0_75, \
                         skd0_77, skd1_75, skd1_77, skf_122, skf_123, \
                         skf_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_12 * sif_123[k]
                   + f_4 * skd0_75[k]
                   - f_5 * skd1_75[k]
                   + f_3 * pc_x[k] * skf_123[k];

        t_184[k] = f_8 * sif_82[k]
                   + f_3 * pc_y[k] * skf_122[k];

        t_185[k] = f_12 * sif_125[k]
                   + f_4 * skd0_77[k]
                   - f_5 * skd1_77[k]
                   + f_3 * pc_x[k] * skf_125[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, pc_x, sif_126, sif_127, sif_128, sif_129, \
                         skf_126, skf_127, skf_128, skf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_12 * sif_126[k]
                   + f_3 * pc_x[k] * skf_126[k];

        t_187[k] = f_12 * sif_127[k]
                   + f_3 * pc_x[k] * skf_127[k];

        t_188[k] = f_12 * sif_128[k]
                   + f_3 * pc_x[k] * skf_128[k];

        t_189[k] = f_12 * sif_129[k]
                   + f_3 * pc_x[k] * skf_129[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, pc_y, pc_z, sif_76, sif_86, sif_88, skd0_75, \
                         skd0_77, skd1_75, skd1_77, skf_126, skf_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = f_8 * sif_86[k]
                   + f_1 * skd0_75[k]
                   - f_2 * skd1_75[k]
                   + f_3 * pc_y[k] * skf_126[k];

        t_191[k] = f_8 * sif_76[k]
                   + f_3 * pc_z[k] * skf_126[k];

        t_192[k] = f_8 * sif_88[k]
                   + f_4 * skd0_77[k]
                   - f_5 * skd1_77[k]
                   + f_3 * pc_y[k] * skf_128[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, t_196, pb_y, pc_y, pc_z, sig0_135, sif_79, \
                         sif_89, sif_90, sig1_135, skd0_77, skd1_77, skf_129, \
                         skf_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_8 * sif_89[k]
                   + f_3 * pc_y[k] * skf_129[k];

        t_194[k] = f_8 * sif_79[k]
                   + f_1 * skd0_77[k]
                   - f_2 * skd1_77[k]
                   + f_3 * pc_z[k] * skf_129[k];

        t_195[k] = pb_y[k] * sig0_135[k]
                   - f_6 * pc_y[k] * sig1_135[k];

        t_196[k] = f_7 * sif_90[k]
                   + f_3 * pc_y[k] * skf_130[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, pb_y, pc_y, pc_z, sig0_138, sig0_140, \
                         sif_80, sif_91, sif_92, sig1_138, sig1_140, skf_130, \
                         skf_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_12 * sif_80[k]
                   + f_3 * pc_z[k] * skf_130[k];

        t_198[k] = pb_y[k] * sig0_138[k]
                   + f_8 * sif_91[k]
                   - f_6 * pc_y[k] * sig1_138[k];

        t_199[k] = f_7 * sif_92[k]
                   + f_3 * pc_y[k] * skf_132[k];

        t_200[k] = pb_y[k] * sig0_140[k]
                   - f_6 * pc_y[k] * sig1_140[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, pc_x, sif_136, sif_137, sif_138, sif_139, \
                         skf_136, skf_137, skf_138, skf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_12 * sif_136[k]
                   + f_3 * pc_x[k] * skf_136[k];

        t_202[k] = f_12 * sif_137[k]
                   + f_3 * pc_x[k] * skf_137[k];

        t_203[k] = f_12 * sif_138[k]
                   + f_3 * pc_x[k] * skf_138[k];

        t_204[k] = f_12 * sif_139[k]
                   + f_3 * pc_x[k] * skf_139[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, pc_y, pc_z, sif_86, sif_96, sif_98, skd0_81, \
                         skd0_83, skd1_81, skd1_83, skf_136, skf_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_7 * sif_96[k]
                   + f_1 * skd0_81[k]
                   - f_2 * skd1_81[k]
                   + f_3 * pc_y[k] * skf_136[k];

        t_206[k] = f_12 * sif_86[k]
                   + f_3 * pc_z[k] * skf_136[k];

        t_207[k] = f_7 * sif_98[k]
                   + f_4 * skd0_83[k]
                   - f_5 * skd1_83[k]
                   + f_3 * pc_y[k] * skf_138[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pb_y, pc_x, pc_y, sig0_149, sif_99, \
                         sif_140, sig1_149, skd0_84, skd1_84, skf_139, \
                         skf_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_7 * sif_99[k]
                   + f_3 * pc_y[k] * skf_139[k];

        t_209[k] = pb_y[k] * sig0_149[k]
                   - f_6 * pc_y[k] * sig1_149[k];

        t_210[k] = f_12 * sif_140[k]
                   + f_1 * skd0_84[k]
                   - f_2 * skd1_84[k]
                   + f_3 * pc_x[k] * skf_140[k];

        t_211[k] = f_3 * pc_y[k] * skf_140[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, pc_x, pc_y, pc_z, sif_90, sif_143, skd0_87, \
                         skd1_87, skf_140, skf_142, skf_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_11 * sif_90[k]
                   + f_3 * pc_z[k] * skf_140[k];

        t_213[k] = f_12 * sif_143[k]
                   + f_4 * skd0_87[k]
                   - f_5 * skd1_87[k]
                   + f_3 * pc_x[k] * skf_143[k];

        t_214[k] = f_3 * pc_y[k] * skf_142[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, t_218, pc_x, sif_145, sif_146, sif_147, sif_148, \
                         skd0_89, skd1_89, skf_145, skf_146, skf_147, \
                         skf_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_12 * sif_145[k]
                   + f_4 * skd0_89[k]
                   - f_5 * skd1_89[k]
                   + f_3 * pc_x[k] * skf_145[k];

        t_216[k] = f_12 * sif_146[k]
                   + f_3 * pc_x[k] * skf_146[k];

        t_217[k] = f_12 * sif_147[k]
                   + f_3 * pc_x[k] * skf_147[k];

        t_218[k] = f_12 * sif_148[k]
                   + f_3 * pc_x[k] * skf_148[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, t_222, pc_x, pc_y, pc_z, sif_96, sif_149, \
                         skd0_87, skd0_89, skd1_87, skd1_89, skf_146, skf_148, \
                         skf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_12 * sif_149[k]
                   + f_3 * pc_x[k] * skf_149[k];

        t_220[k] = f_1 * skd0_87[k]
                   - f_2 * skd1_87[k]
                   + f_3 * pc_y[k] * skf_146[k];

        t_221[k] = f_11 * sif_96[k]
                   + f_3 * pc_z[k] * skf_146[k];

        t_222[k] = f_4 * skd0_89[k]
                   - f_5 * skd1_89[k]
                   + f_3 * pc_y[k] * skf_148[k];
    }

#pragma omp simd aligned(t_223, t_224, t_225, t_226, pc_x, pc_y, pc_z, sif_99, sif_100, \
                         sif_150, skd0_89, skd0_90, skd1_89, skd1_90, skf_149, \
                         skf_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_223[k] = f_3 * pc_y[k] * skf_149[k];

        t_224[k] = f_11 * sif_99[k]
                   + f_1 * skd0_89[k]
                   - f_2 * skd1_89[k]
                   + f_3 * pc_z[k] * skf_149[k];

        t_225[k] = f_8 * sif_150[k]
                   + f_1 * skd0_90[k]
                   - f_2 * skd1_90[k]
                   + f_3 * pc_x[k] * skf_150[k];

        t_226[k] = f_10 * sif_100[k]
                   + f_3 * pc_y[k] * skf_150[k];
    }

#pragma omp simd aligned(t_227, t_228, t_229, pc_x, pc_y, pc_z, sif_102, sif_153, skd0_93, \
                         skd1_93, skf_150, skf_152, skf_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_227[k] = f_3 * pc_z[k] * skf_150[k];

        t_228[k] = f_8 * sif_153[k]
                   + f_4 * skd0_93[k]
                   - f_5 * skd1_93[k]
                   + f_3 * pc_x[k] * skf_153[k];

        t_229[k] = f_10 * sif_102[k]
                   + f_3 * pc_y[k] * skf_152[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, t_233, pc_x, sif_155, sif_156, sif_157, sif_158, \
                         skd0_95, skd1_95, skf_155, skf_156, skf_157, \
                         skf_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = f_8 * sif_155[k]
                   + f_4 * skd0_95[k]
                   - f_5 * skd1_95[k]
                   + f_3 * pc_x[k] * skf_155[k];

        t_231[k] = f_8 * sif_156[k]
                   + f_3 * pc_x[k] * skf_156[k];

        t_232[k] = f_8 * sif_157[k]
                   + f_3 * pc_x[k] * skf_157[k];

        t_233[k] = f_8 * sif_158[k]
                   + f_3 * pc_x[k] * skf_158[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, pc_x, pc_y, pc_z, sif_106, sif_159, skd0_93, \
                         skd1_93, skf_156, skf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_8 * sif_159[k]
                   + f_3 * pc_x[k] * skf_159[k];

        t_235[k] = f_10 * sif_106[k]
                   + f_1 * skd0_93[k]
                   - f_2 * skd1_93[k]
                   + f_3 * pc_y[k] * skf_156[k];

        t_236[k] = f_3 * pc_z[k] * skf_156[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, t_240, pb_z, pc_y, pc_z, sig0_150, sif_108, \
                         sif_109, sig1_150, skd0_95, skd1_95, skf_158, \
                         skf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_10 * sif_108[k]
                   + f_4 * skd0_95[k]
                   - f_5 * skd1_95[k]
                   + f_3 * pc_y[k] * skf_158[k];

        t_238[k] = f_10 * sif_109[k]
                   + f_3 * pc_y[k] * skf_159[k];

        t_239[k] = f_1 * skd0_95[k]
                   - f_2 * skd1_95[k]
                   + f_3 * pc_z[k] * skf_159[k];

        t_240[k] = pb_z[k] * sig0_150[k]
                   - f_6 * pc_z[k] * sig1_150[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, pb_z, pc_y, pc_z, sig0_153, sif_100, \
                         sif_110, sif_112, sig1_153, skf_160, skf_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = f_11 * sif_110[k]
                   + f_3 * pc_y[k] * skf_160[k];

        t_242[k] = f_7 * sif_100[k]
                   + f_3 * pc_z[k] * skf_160[k];

        t_243[k] = pb_z[k] * sig0_153[k]
                   - f_6 * pc_z[k] * sig1_153[k];

        t_244[k] = f_11 * sif_112[k]
                   + f_3 * pc_y[k] * skf_162[k];
    }
}

static auto
compute_prim_skg_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sig0,
                                                          const size_t sif, const size_t sig1,
                                                          const size_t skd0, const size_t skd1,
                                                          const size_t skf, const size_t ncols,
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
    const auto f_9 = 3.0 / q;
    const auto f_10 = 2.5 / q;
    const auto f_11 = 2.0 / q;
    const auto f_12 = 1.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sig0_160 = buffer.data(sig0 + 160);
    const auto *sig0_210 = buffer.data(sig0 + 210);
    const auto *sig0_213 = buffer.data(sig0 + 213);
    const auto *sig0_215 = buffer.data(sig0 + 215);
    const auto *sig0_224 = buffer.data(sig0 + 224);
    const auto *sig0_225 = buffer.data(sig0 + 225);
    const auto *sig0_228 = buffer.data(sig0 + 228);
    const auto *sig0_315 = buffer.data(sig0 + 315);
    const auto *sig0_318 = buffer.data(sig0 + 318);
    const auto *sig0_320 = buffer.data(sig0 + 320);
    const auto *sig0_325 = buffer.data(sig0 + 325);
    const auto *sig0_327 = buffer.data(sig0 + 327);
    const auto *sig0_329 = buffer.data(sig0 + 329);
    const auto *sig0_335 = buffer.data(sig0 + 335);
    const auto *sig0_340 = buffer.data(sig0 + 340);
    const auto *sig0_342 = buffer.data(sig0 + 342);
    const auto *sig0_344 = buffer.data(sig0 + 344);
    const auto *sig0_345 = buffer.data(sig0 + 345);
    const auto *sig0_348 = buffer.data(sig0 + 348);
    const auto *sig0_350 = buffer.data(sig0 + 350);
    const auto *sig0_355 = buffer.data(sig0 + 355);
    const auto *sig0_357 = buffer.data(sig0 + 357);
    const auto *sig0_359 = buffer.data(sig0 + 359);
    const auto *sig0_360 = buffer.data(sig0 + 360);
    const auto *sig0_363 = buffer.data(sig0 + 363);

    const auto *sif_106 = buffer.data(sif + 106);
    const auto *sif_109 = buffer.data(sif + 109);
    const auto *sif_110 = buffer.data(sif + 110);
    const auto *sif_116 = buffer.data(sif + 116);
    const auto *sif_118 = buffer.data(sif + 118);
    const auto *sif_119 = buffer.data(sif + 119);
    const auto *sif_120 = buffer.data(sif + 120);
    const auto *sif_122 = buffer.data(sif + 122);
    const auto *sif_126 = buffer.data(sif + 126);
    const auto *sif_128 = buffer.data(sif + 128);
    const auto *sif_129 = buffer.data(sif + 129);
    const auto *sif_130 = buffer.data(sif + 130);
    const auto *sif_132 = buffer.data(sif + 132);
    const auto *sif_136 = buffer.data(sif + 136);
    const auto *sif_138 = buffer.data(sif + 138);
    const auto *sif_139 = buffer.data(sif + 139);
    const auto *sif_140 = buffer.data(sif + 140);
    const auto *sif_141 = buffer.data(sif + 141);
    const auto *sif_142 = buffer.data(sif + 142);
    const auto *sif_146 = buffer.data(sif + 146);
    const auto *sif_148 = buffer.data(sif + 148);
    const auto *sif_149 = buffer.data(sif + 149);
    const auto *sif_150 = buffer.data(sif + 150);
    const auto *sif_152 = buffer.data(sif + 152);
    const auto *sif_156 = buffer.data(sif + 156);
    const auto *sif_159 = buffer.data(sif + 159);
    const auto *sif_160 = buffer.data(sif + 160);
    const auto *sif_162 = buffer.data(sif + 162);
    const auto *sif_165 = buffer.data(sif + 165);
    const auto *sif_166 = buffer.data(sif + 166);
    const auto *sif_167 = buffer.data(sif + 167);
    const auto *sif_168 = buffer.data(sif + 168);
    const auto *sif_169 = buffer.data(sif + 169);
    const auto *sif_170 = buffer.data(sif + 170);
    const auto *sif_172 = buffer.data(sif + 172);
    const auto *sif_173 = buffer.data(sif + 173);
    const auto *sif_175 = buffer.data(sif + 175);
    const auto *sif_176 = buffer.data(sif + 176);
    const auto *sif_177 = buffer.data(sif + 177);
    const auto *sif_178 = buffer.data(sif + 178);
    const auto *sif_179 = buffer.data(sif + 179);
    const auto *sif_180 = buffer.data(sif + 180);
    const auto *sif_182 = buffer.data(sif + 182);
    const auto *sif_183 = buffer.data(sif + 183);
    const auto *sif_185 = buffer.data(sif + 185);
    const auto *sif_186 = buffer.data(sif + 186);
    const auto *sif_187 = buffer.data(sif + 187);
    const auto *sif_188 = buffer.data(sif + 188);
    const auto *sif_189 = buffer.data(sif + 189);
    const auto *sif_196 = buffer.data(sif + 196);
    const auto *sif_197 = buffer.data(sif + 197);
    const auto *sif_198 = buffer.data(sif + 198);
    const auto *sif_199 = buffer.data(sif + 199);
    const auto *sif_200 = buffer.data(sif + 200);
    const auto *sif_203 = buffer.data(sif + 203);
    const auto *sif_205 = buffer.data(sif + 205);
    const auto *sif_206 = buffer.data(sif + 206);
    const auto *sif_207 = buffer.data(sif + 207);
    const auto *sif_208 = buffer.data(sif + 208);
    const auto *sif_209 = buffer.data(sif + 209);
    const auto *sif_210 = buffer.data(sif + 210);
    const auto *sif_213 = buffer.data(sif + 213);
    const auto *sif_215 = buffer.data(sif + 215);
    const auto *sif_216 = buffer.data(sif + 216);
    const auto *sif_217 = buffer.data(sif + 217);
    const auto *sif_218 = buffer.data(sif + 218);
    const auto *sif_219 = buffer.data(sif + 219);
    const auto *sif_225 = buffer.data(sif + 225);
    const auto *sif_226 = buffer.data(sif + 226);
    const auto *sif_227 = buffer.data(sif + 227);
    const auto *sif_228 = buffer.data(sif + 228);
    const auto *sif_229 = buffer.data(sif + 229);
    const auto *sif_230 = buffer.data(sif + 230);
    const auto *sif_233 = buffer.data(sif + 233);
    const auto *sif_235 = buffer.data(sif + 235);
    const auto *sif_236 = buffer.data(sif + 236);
    const auto *sif_237 = buffer.data(sif + 237);
    const auto *sif_238 = buffer.data(sif + 238);
    const auto *sif_239 = buffer.data(sif + 239);
    const auto *sif_240 = buffer.data(sif + 240);
    const auto *sif_243 = buffer.data(sif + 243);

    const auto *sig1_160 = buffer.data(sig1 + 160);
    const auto *sig1_210 = buffer.data(sig1 + 210);
    const auto *sig1_213 = buffer.data(sig1 + 213);
    const auto *sig1_215 = buffer.data(sig1 + 215);
    const auto *sig1_224 = buffer.data(sig1 + 224);
    const auto *sig1_225 = buffer.data(sig1 + 225);
    const auto *sig1_228 = buffer.data(sig1 + 228);
    const auto *sig1_315 = buffer.data(sig1 + 315);
    const auto *sig1_318 = buffer.data(sig1 + 318);
    const auto *sig1_320 = buffer.data(sig1 + 320);
    const auto *sig1_325 = buffer.data(sig1 + 325);
    const auto *sig1_327 = buffer.data(sig1 + 327);
    const auto *sig1_329 = buffer.data(sig1 + 329);
    const auto *sig1_335 = buffer.data(sig1 + 335);
    const auto *sig1_340 = buffer.data(sig1 + 340);
    const auto *sig1_342 = buffer.data(sig1 + 342);
    const auto *sig1_344 = buffer.data(sig1 + 344);
    const auto *sig1_345 = buffer.data(sig1 + 345);
    const auto *sig1_348 = buffer.data(sig1 + 348);
    const auto *sig1_350 = buffer.data(sig1 + 350);
    const auto *sig1_355 = buffer.data(sig1 + 355);
    const auto *sig1_357 = buffer.data(sig1 + 357);
    const auto *sig1_359 = buffer.data(sig1 + 359);
    const auto *sig1_360 = buffer.data(sig1 + 360);
    const auto *sig1_363 = buffer.data(sig1 + 363);

    const auto *skd0_101 = buffer.data(skd0 + 101);
    const auto *skd0_102 = buffer.data(skd0 + 102);
    const auto *skd0_105 = buffer.data(skd0 + 105);
    const auto *skd0_107 = buffer.data(skd0 + 107);
    const auto *skd0_108 = buffer.data(skd0 + 108);
    const auto *skd0_111 = buffer.data(skd0 + 111);
    const auto *skd0_113 = buffer.data(skd0 + 113);
    const auto *skd0_117 = buffer.data(skd0 + 117);
    const auto *skd0_119 = buffer.data(skd0 + 119);
    const auto *skd0_120 = buffer.data(skd0 + 120);
    const auto *skd0_123 = buffer.data(skd0 + 123);
    const auto *skd0_125 = buffer.data(skd0 + 125);

    const auto *skd1_101 = buffer.data(skd1 + 101);
    const auto *skd1_102 = buffer.data(skd1 + 102);
    const auto *skd1_105 = buffer.data(skd1 + 105);
    const auto *skd1_107 = buffer.data(skd1 + 107);
    const auto *skd1_108 = buffer.data(skd1 + 108);
    const auto *skd1_111 = buffer.data(skd1 + 111);
    const auto *skd1_113 = buffer.data(skd1 + 113);
    const auto *skd1_117 = buffer.data(skd1 + 117);
    const auto *skd1_119 = buffer.data(skd1 + 119);
    const auto *skd1_120 = buffer.data(skd1 + 120);
    const auto *skd1_123 = buffer.data(skd1 + 123);
    const auto *skd1_125 = buffer.data(skd1 + 125);

    const auto *skf_165 = buffer.data(skf + 165);
    const auto *skf_166 = buffer.data(skf + 166);
    const auto *skf_167 = buffer.data(skf + 167);
    const auto *skf_168 = buffer.data(skf + 168);
    const auto *skf_169 = buffer.data(skf + 169);
    const auto *skf_170 = buffer.data(skf + 170);
    const auto *skf_172 = buffer.data(skf + 172);
    const auto *skf_173 = buffer.data(skf + 173);
    const auto *skf_175 = buffer.data(skf + 175);
    const auto *skf_176 = buffer.data(skf + 176);
    const auto *skf_177 = buffer.data(skf + 177);
    const auto *skf_178 = buffer.data(skf + 178);
    const auto *skf_179 = buffer.data(skf + 179);
    const auto *skf_180 = buffer.data(skf + 180);
    const auto *skf_182 = buffer.data(skf + 182);
    const auto *skf_183 = buffer.data(skf + 183);
    const auto *skf_185 = buffer.data(skf + 185);
    const auto *skf_186 = buffer.data(skf + 186);
    const auto *skf_187 = buffer.data(skf + 187);
    const auto *skf_188 = buffer.data(skf + 188);
    const auto *skf_189 = buffer.data(skf + 189);
    const auto *skf_190 = buffer.data(skf + 190);
    const auto *skf_192 = buffer.data(skf + 192);
    const auto *skf_196 = buffer.data(skf + 196);
    const auto *skf_197 = buffer.data(skf + 197);
    const auto *skf_198 = buffer.data(skf + 198);
    const auto *skf_199 = buffer.data(skf + 199);
    const auto *skf_200 = buffer.data(skf + 200);
    const auto *skf_202 = buffer.data(skf + 202);
    const auto *skf_203 = buffer.data(skf + 203);
    const auto *skf_205 = buffer.data(skf + 205);
    const auto *skf_206 = buffer.data(skf + 206);
    const auto *skf_207 = buffer.data(skf + 207);
    const auto *skf_208 = buffer.data(skf + 208);
    const auto *skf_209 = buffer.data(skf + 209);
    const auto *skf_210 = buffer.data(skf + 210);
    const auto *skf_212 = buffer.data(skf + 212);
    const auto *skf_216 = buffer.data(skf + 216);
    const auto *skf_217 = buffer.data(skf + 217);
    const auto *skf_218 = buffer.data(skf + 218);
    const auto *skf_219 = buffer.data(skf + 219);
    const auto *skf_220 = buffer.data(skf + 220);
    const auto *skf_222 = buffer.data(skf + 222);
    const auto *skf_226 = buffer.data(skf + 226);
    const auto *skf_227 = buffer.data(skf + 227);
    const auto *skf_228 = buffer.data(skf + 228);
    const auto *skf_229 = buffer.data(skf + 229);
    const auto *skf_230 = buffer.data(skf + 230);
    const auto *skf_232 = buffer.data(skf + 232);
    const auto *skf_236 = buffer.data(skf + 236);
    const auto *skf_237 = buffer.data(skf + 237);
    const auto *skf_238 = buffer.data(skf + 238);
    const auto *skf_239 = buffer.data(skf + 239);
    const auto *skf_240 = buffer.data(skf + 240);
    const auto *skf_242 = buffer.data(skf + 242);

#pragma omp simd aligned(t_245, t_246, t_247, t_248, pc_x, sif_165, sif_166, sif_167, sif_168, \
                         skd0_101, skd1_101, skf_165, skf_166, skf_167, \
                         skf_168 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_8 * sif_165[k]
                   + f_4 * skd0_101[k]
                   - f_5 * skd1_101[k]
                   + f_3 * pc_x[k] * skf_165[k];

        t_246[k] = f_8 * sif_166[k]
                   + f_3 * pc_x[k] * skf_166[k];

        t_247[k] = f_8 * sif_167[k]
                   + f_3 * pc_x[k] * skf_167[k];

        t_248[k] = f_8 * sif_168[k]
                   + f_3 * pc_x[k] * skf_168[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, pb_z, pc_x, pc_z, sig0_160, sif_106, sif_169, \
                         sig1_160, skf_166, skf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_8 * sif_169[k]
                   + f_3 * pc_x[k] * skf_169[k];

        t_250[k] = pb_z[k] * sig0_160[k]
                   - f_6 * pc_z[k] * sig1_160[k];

        t_251[k] = f_7 * sif_106[k]
                   + f_3 * pc_z[k] * skf_166[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, pc_y, pc_z, sif_109, sif_118, sif_119, skd0_101, \
                         skd1_101, skf_168, skf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_11 * sif_118[k]
                   + f_4 * skd0_101[k]
                   - f_5 * skd1_101[k]
                   + f_3 * pc_y[k] * skf_168[k];

        t_253[k] = f_11 * sif_119[k]
                   + f_3 * pc_y[k] * skf_169[k];

        t_254[k] = f_7 * sif_109[k]
                   + f_1 * skd0_101[k]
                   - f_2 * skd1_101[k]
                   + f_3 * pc_z[k] * skf_169[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, pc_x, pc_y, pc_z, sif_110, sif_120, sif_170, \
                         skd0_102, skd1_102, skf_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_8 * sif_170[k]
                   + f_1 * skd0_102[k]
                   - f_2 * skd1_102[k]
                   + f_3 * pc_x[k] * skf_170[k];

        t_256[k] = f_12 * sif_120[k]
                   + f_3 * pc_y[k] * skf_170[k];

        t_257[k] = f_8 * sif_110[k]
                   + f_3 * pc_z[k] * skf_170[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, pc_x, pc_y, sif_122, sif_173, sif_175, skd0_105, \
                         skd0_107, skd1_105, skd1_107, skf_172, skf_173, \
                         skf_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_8 * sif_173[k]
                   + f_4 * skd0_105[k]
                   - f_5 * skd1_105[k]
                   + f_3 * pc_x[k] * skf_173[k];

        t_259[k] = f_12 * sif_122[k]
                   + f_3 * pc_y[k] * skf_172[k];

        t_260[k] = f_8 * sif_175[k]
                   + f_4 * skd0_107[k]
                   - f_5 * skd1_107[k]
                   + f_3 * pc_x[k] * skf_175[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, pc_x, sif_176, sif_177, sif_178, sif_179, \
                         skf_176, skf_177, skf_178, skf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_8 * sif_176[k]
                   + f_3 * pc_x[k] * skf_176[k];

        t_262[k] = f_8 * sif_177[k]
                   + f_3 * pc_x[k] * skf_177[k];

        t_263[k] = f_8 * sif_178[k]
                   + f_3 * pc_x[k] * skf_178[k];

        t_264[k] = f_8 * sif_179[k]
                   + f_3 * pc_x[k] * skf_179[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, pc_y, pc_z, sif_116, sif_126, sif_128, skd0_105, \
                         skd0_107, skd1_105, skd1_107, skf_176, \
                         skf_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_12 * sif_126[k]
                   + f_1 * skd0_105[k]
                   - f_2 * skd1_105[k]
                   + f_3 * pc_y[k] * skf_176[k];

        t_266[k] = f_8 * sif_116[k]
                   + f_3 * pc_z[k] * skf_176[k];

        t_267[k] = f_12 * sif_128[k]
                   + f_4 * skd0_107[k]
                   - f_5 * skd1_107[k]
                   + f_3 * pc_y[k] * skf_178[k];
    }

#pragma omp simd aligned(t_268, t_269, t_270, pc_x, pc_y, pc_z, sif_119, sif_129, sif_180, \
                         skd0_107, skd0_108, skd1_107, skd1_108, skf_179, \
                         skf_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_12 * sif_129[k]
                   + f_3 * pc_y[k] * skf_179[k];

        t_269[k] = f_8 * sif_119[k]
                   + f_1 * skd0_107[k]
                   - f_2 * skd1_107[k]
                   + f_3 * pc_z[k] * skf_179[k];

        t_270[k] = f_8 * sif_180[k]
                   + f_1 * skd0_108[k]
                   - f_2 * skd1_108[k]
                   + f_3 * pc_x[k] * skf_180[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, t_274, pc_x, pc_y, pc_z, sif_120, sif_130, \
                         sif_132, sif_183, skd0_111, skd1_111, skf_180, skf_182, \
                         skf_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_8 * sif_130[k]
                   + f_3 * pc_y[k] * skf_180[k];

        t_272[k] = f_12 * sif_120[k]
                   + f_3 * pc_z[k] * skf_180[k];

        t_273[k] = f_8 * sif_183[k]
                   + f_4 * skd0_111[k]
                   - f_5 * skd1_111[k]
                   + f_3 * pc_x[k] * skf_183[k];

        t_274[k] = f_8 * sif_132[k]
                   + f_3 * pc_y[k] * skf_182[k];
    }

#pragma omp simd aligned(t_275, t_276, t_277, t_278, pc_x, sif_185, sif_186, sif_187, sif_188, \
                         skd0_113, skd1_113, skf_185, skf_186, skf_187, \
                         skf_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_275[k] = f_8 * sif_185[k]
                   + f_4 * skd0_113[k]
                   - f_5 * skd1_113[k]
                   + f_3 * pc_x[k] * skf_185[k];

        t_276[k] = f_8 * sif_186[k]
                   + f_3 * pc_x[k] * skf_186[k];

        t_277[k] = f_8 * sif_187[k]
                   + f_3 * pc_x[k] * skf_187[k];

        t_278[k] = f_8 * sif_188[k]
                   + f_3 * pc_x[k] * skf_188[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, pc_x, pc_y, pc_z, sif_126, sif_136, sif_189, \
                         skd0_111, skd1_111, skf_186, skf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = f_8 * sif_189[k]
                   + f_3 * pc_x[k] * skf_189[k];

        t_280[k] = f_8 * sif_136[k]
                   + f_1 * skd0_111[k]
                   - f_2 * skd1_111[k]
                   + f_3 * pc_y[k] * skf_186[k];

        t_281[k] = f_12 * sif_126[k]
                   + f_3 * pc_z[k] * skf_186[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, t_285, pb_y, pc_y, pc_z, sig0_210, sif_129, \
                         sif_138, sif_139, sig1_210, skd0_113, skd1_113, skf_188, \
                         skf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = f_8 * sif_138[k]
                   + f_4 * skd0_113[k]
                   - f_5 * skd1_113[k]
                   + f_3 * pc_y[k] * skf_188[k];

        t_283[k] = f_8 * sif_139[k]
                   + f_3 * pc_y[k] * skf_189[k];

        t_284[k] = f_12 * sif_129[k]
                   + f_1 * skd0_113[k]
                   - f_2 * skd1_113[k]
                   + f_3 * pc_z[k] * skf_189[k];

        t_285[k] = pb_y[k] * sig0_210[k]
                   - f_6 * pc_y[k] * sig1_210[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pb_y, pc_y, pc_z, sig0_213, sif_130, \
                         sif_140, sif_141, sif_142, sig1_213, skf_190, \
                         skf_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_7 * sif_140[k]
                   + f_3 * pc_y[k] * skf_190[k];

        t_287[k] = f_11 * sif_130[k]
                   + f_3 * pc_z[k] * skf_190[k];

        t_288[k] = pb_y[k] * sig0_213[k]
                   + f_8 * sif_141[k]
                   - f_6 * pc_y[k] * sig1_213[k];

        t_289[k] = f_7 * sif_142[k]
                   + f_3 * pc_y[k] * skf_192[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, pb_y, pc_x, pc_y, sig0_215, sif_196, \
                         sif_197, sif_198, sig1_215, skf_196, skf_197, \
                         skf_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = pb_y[k] * sig0_215[k]
                   - f_6 * pc_y[k] * sig1_215[k];

        t_291[k] = f_8 * sif_196[k]
                   + f_3 * pc_x[k] * skf_196[k];

        t_292[k] = f_8 * sif_197[k]
                   + f_3 * pc_x[k] * skf_197[k];

        t_293[k] = f_8 * sif_198[k]
                   + f_3 * pc_x[k] * skf_198[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, pc_x, pc_y, pc_z, sif_136, sif_146, sif_199, \
                         skd0_117, skd1_117, skf_196, skf_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_8 * sif_199[k]
                   + f_3 * pc_x[k] * skf_199[k];

        t_295[k] = f_7 * sif_146[k]
                   + f_1 * skd0_117[k]
                   - f_2 * skd1_117[k]
                   + f_3 * pc_y[k] * skf_196[k];

        t_296[k] = f_11 * sif_136[k]
                   + f_3 * pc_z[k] * skf_196[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, pb_y, pc_y, sig0_224, sif_148, sif_149, \
                         sig1_224, skd0_119, skd1_119, skf_198, \
                         skf_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_7 * sif_148[k]
                   + f_4 * skd0_119[k]
                   - f_5 * skd1_119[k]
                   + f_3 * pc_y[k] * skf_198[k];

        t_298[k] = f_7 * sif_149[k]
                   + f_3 * pc_y[k] * skf_199[k];

        t_299[k] = pb_y[k] * sig0_224[k]
                   - f_6 * pc_y[k] * sig1_224[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, pc_x, pc_y, pc_z, sif_140, sif_200, \
                         sif_203, skd0_120, skd0_123, skd1_120, skd1_123, skf_200, \
                         skf_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_8 * sif_200[k]
                   + f_1 * skd0_120[k]
                   - f_2 * skd1_120[k]
                   + f_3 * pc_x[k] * skf_200[k];

        t_301[k] = f_3 * pc_y[k] * skf_200[k];

        t_302[k] = f_10 * sif_140[k]
                   + f_3 * pc_z[k] * skf_200[k];

        t_303[k] = f_8 * sif_203[k]
                   + f_4 * skd0_123[k]
                   - f_5 * skd1_123[k]
                   + f_3 * pc_x[k] * skf_203[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, pc_x, pc_y, sif_205, sif_206, sif_207, \
                         skd0_125, skd1_125, skf_202, skf_205, skf_206, \
                         skf_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = f_3 * pc_y[k] * skf_202[k];

        t_305[k] = f_8 * sif_205[k]
                   + f_4 * skd0_125[k]
                   - f_5 * skd1_125[k]
                   + f_3 * pc_x[k] * skf_205[k];

        t_306[k] = f_8 * sif_206[k]
                   + f_3 * pc_x[k] * skf_206[k];

        t_307[k] = f_8 * sif_207[k]
                   + f_3 * pc_x[k] * skf_207[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, pc_x, pc_y, pc_z, sif_146, sif_208, \
                         sif_209, skd0_123, skd1_123, skf_206, skf_208, \
                         skf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = f_8 * sif_208[k]
                   + f_3 * pc_x[k] * skf_208[k];

        t_309[k] = f_8 * sif_209[k]
                   + f_3 * pc_x[k] * skf_209[k];

        t_310[k] = f_1 * skd0_123[k]
                   - f_2 * skd1_123[k]
                   + f_3 * pc_y[k] * skf_206[k];

        t_311[k] = f_10 * sif_146[k]
                   + f_3 * pc_z[k] * skf_206[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, t_315, pb_x, pc_x, pc_y, pc_z, sig0_315, \
                         sif_149, sif_210, sig1_315, skd0_125, skd1_125, skf_208, \
                         skf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_4 * skd0_125[k]
                   - f_5 * skd1_125[k]
                   + f_3 * pc_y[k] * skf_208[k];

        t_313[k] = f_3 * pc_y[k] * skf_209[k];

        t_314[k] = f_10 * sif_149[k]
                   + f_1 * skd0_125[k]
                   - f_2 * skd1_125[k]
                   + f_3 * pc_z[k] * skf_209[k];

        t_315[k] = pb_x[k] * sig0_315[k]
                   + f_11 * sif_210[k]
                   - f_6 * pc_x[k] * sig1_315[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, t_319, pb_x, pc_x, pc_y, pc_z, sig0_318, \
                         sif_150, sif_152, sif_213, sig1_318, skf_210, \
                         skf_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = f_9 * sif_150[k]
                   + f_3 * pc_y[k] * skf_210[k];

        t_317[k] = f_3 * pc_z[k] * skf_210[k];

        t_318[k] = pb_x[k] * sig0_318[k]
                   + f_8 * sif_213[k]
                   - f_6 * pc_x[k] * sig1_318[k];

        t_319[k] = f_9 * sif_152[k]
                   + f_3 * pc_y[k] * skf_212[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, pb_x, pc_x, sig0_320, sif_215, sif_216, \
                         sif_217, sif_218, sig1_320, skf_216, skf_217, \
                         skf_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = pb_x[k] * sig0_320[k]
                   + f_8 * sif_215[k]
                   - f_6 * pc_x[k] * sig1_320[k];

        t_321[k] = f_7 * sif_216[k]
                   + f_3 * pc_x[k] * skf_216[k];

        t_322[k] = f_7 * sif_217[k]
                   + f_3 * pc_x[k] * skf_217[k];

        t_323[k] = f_7 * sif_218[k]
                   + f_3 * pc_x[k] * skf_218[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, pb_x, pc_x, pc_z, sig0_325, sig0_327, \
                         sif_219, sig1_325, sig1_327, skf_216, \
                         skf_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_7 * sif_219[k]
                   + f_3 * pc_x[k] * skf_219[k];

        t_325[k] = pb_x[k] * sig0_325[k]
                   - f_6 * pc_x[k] * sig1_325[k];

        t_326[k] = f_3 * pc_z[k] * skf_216[k];

        t_327[k] = pb_x[k] * sig0_327[k]
                   - f_6 * pc_x[k] * sig1_327[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, pb_x, pb_z, pc_x, pc_y, pc_z, sig0_225, \
                         sig0_329, sif_159, sig1_225, sig1_329, \
                         skf_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = f_9 * sif_159[k]
                   + f_3 * pc_y[k] * skf_219[k];

        t_329[k] = pb_x[k] * sig0_329[k]
                   - f_6 * pc_x[k] * sig1_329[k];

        t_330[k] = pb_z[k] * sig0_225[k]
                   - f_6 * pc_z[k] * sig1_225[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, t_334, pb_z, pc_y, pc_z, sig0_228, sif_150, \
                         sif_160, sif_162, sig1_228, skf_220, skf_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = f_10 * sif_160[k]
                   + f_3 * pc_y[k] * skf_220[k];

        t_332[k] = f_7 * sif_150[k]
                   + f_3 * pc_z[k] * skf_220[k];

        t_333[k] = pb_z[k] * sig0_228[k]
                   - f_6 * pc_z[k] * sig1_228[k];

        t_334[k] = f_10 * sif_162[k]
                   + f_3 * pc_y[k] * skf_222[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, pb_x, pc_x, sig0_335, sif_225, sif_226, \
                         sif_227, sif_228, sig1_335, skf_226, skf_227, \
                         skf_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = pb_x[k] * sig0_335[k]
                   + f_8 * sif_225[k]
                   - f_6 * pc_x[k] * sig1_335[k];

        t_336[k] = f_7 * sif_226[k]
                   + f_3 * pc_x[k] * skf_226[k];

        t_337[k] = f_7 * sif_227[k]
                   + f_3 * pc_x[k] * skf_227[k];

        t_338[k] = f_7 * sif_228[k]
                   + f_3 * pc_x[k] * skf_228[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, t_342, pb_x, pc_x, pc_z, sig0_340, sig0_342, \
                         sif_156, sif_229, sig1_340, sig1_342, skf_226, \
                         skf_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_7 * sif_229[k]
                   + f_3 * pc_x[k] * skf_229[k];

        t_340[k] = pb_x[k] * sig0_340[k]
                   - f_6 * pc_x[k] * sig1_340[k];

        t_341[k] = f_7 * sif_156[k]
                   + f_3 * pc_z[k] * skf_226[k];

        t_342[k] = pb_x[k] * sig0_342[k]
                   - f_6 * pc_x[k] * sig1_342[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, t_346, pb_x, pc_x, pc_y, sig0_344, sig0_345, \
                         sif_169, sif_170, sif_230, sig1_344, sig1_345, skf_229, \
                         skf_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = f_10 * sif_169[k]
                   + f_3 * pc_y[k] * skf_229[k];

        t_344[k] = pb_x[k] * sig0_344[k]
                   - f_6 * pc_x[k] * sig1_344[k];

        t_345[k] = pb_x[k] * sig0_345[k]
                   + f_11 * sif_230[k]
                   - f_6 * pc_x[k] * sig1_345[k];

        t_346[k] = f_11 * sif_170[k]
                   + f_3 * pc_y[k] * skf_230[k];
    }

#pragma omp simd aligned(t_347, t_348, t_349, pb_x, pc_x, pc_y, pc_z, sig0_348, sif_160, \
                         sif_172, sif_233, sig1_348, skf_230, skf_232 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_347[k] = f_8 * sif_160[k]
                   + f_3 * pc_z[k] * skf_230[k];

        t_348[k] = pb_x[k] * sig0_348[k]
                   + f_8 * sif_233[k]
                   - f_6 * pc_x[k] * sig1_348[k];

        t_349[k] = f_11 * sif_172[k]
                   + f_3 * pc_y[k] * skf_232[k];
    }

#pragma omp simd aligned(t_350, t_351, t_352, t_353, pb_x, pc_x, sig0_350, sif_235, sif_236, \
                         sif_237, sif_238, sig1_350, skf_236, skf_237, \
                         skf_238 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_350[k] = pb_x[k] * sig0_350[k]
                   + f_8 * sif_235[k]
                   - f_6 * pc_x[k] * sig1_350[k];

        t_351[k] = f_7 * sif_236[k]
                   + f_3 * pc_x[k] * skf_236[k];

        t_352[k] = f_7 * sif_237[k]
                   + f_3 * pc_x[k] * skf_237[k];

        t_353[k] = f_7 * sif_238[k]
                   + f_3 * pc_x[k] * skf_238[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, t_357, pb_x, pc_x, pc_z, sig0_355, sig0_357, \
                         sif_166, sif_239, sig1_355, sig1_357, skf_236, \
                         skf_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_7 * sif_239[k]
                   + f_3 * pc_x[k] * skf_239[k];

        t_355[k] = pb_x[k] * sig0_355[k]
                   - f_6 * pc_x[k] * sig1_355[k];

        t_356[k] = f_8 * sif_166[k]
                   + f_3 * pc_z[k] * skf_236[k];

        t_357[k] = pb_x[k] * sig0_357[k]
                   - f_6 * pc_x[k] * sig1_357[k];
    }

#pragma omp simd aligned(t_358, t_359, t_360, t_361, pb_x, pc_x, pc_y, sig0_359, sig0_360, \
                         sif_179, sif_180, sif_240, sig1_359, sig1_360, skf_239, \
                         skf_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_358[k] = f_11 * sif_179[k]
                   + f_3 * pc_y[k] * skf_239[k];

        t_359[k] = pb_x[k] * sig0_359[k]
                   - f_6 * pc_x[k] * sig1_359[k];

        t_360[k] = pb_x[k] * sig0_360[k]
                   + f_11 * sif_240[k]
                   - f_6 * pc_x[k] * sig1_360[k];

        t_361[k] = f_12 * sif_180[k]
                   + f_3 * pc_y[k] * skf_240[k];
    }

#pragma omp simd aligned(t_362, t_363, t_364, pb_x, pc_x, pc_y, pc_z, sig0_363, sif_170, \
                         sif_182, sif_243, sig1_363, skf_240, skf_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_362[k] = f_12 * sif_170[k]
                   + f_3 * pc_z[k] * skf_240[k];

        t_363[k] = pb_x[k] * sig0_363[k]
                   + f_8 * sif_243[k]
                   - f_6 * pc_x[k] * sig1_363[k];

        t_364[k] = f_12 * sif_182[k]
                   + f_3 * pc_y[k] * skf_242[k];
    }
}

static auto
compute_prim_skg_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sig0,
                                                          const size_t sif, const size_t sig1,
                                                          const size_t skd0, const size_t skd1,
                                                          const size_t skf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.0 / q;
    const auto f_10 = 2.5 / q;
    const auto f_11 = 2.0 / q;
    const auto f_12 = 1.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sig0_300 = buffer.data(sig0 + 300);
    const auto *sig0_305 = buffer.data(sig0 + 305);
    const auto *sig0_315 = buffer.data(sig0 + 315);
    const auto *sig0_318 = buffer.data(sig0 + 318);
    const auto *sig0_325 = buffer.data(sig0 + 325);
    const auto *sig0_327 = buffer.data(sig0 + 327);
    const auto *sig0_365 = buffer.data(sig0 + 365);
    const auto *sig0_370 = buffer.data(sig0 + 370);
    const auto *sig0_372 = buffer.data(sig0 + 372);
    const auto *sig0_374 = buffer.data(sig0 + 374);
    const auto *sig0_375 = buffer.data(sig0 + 375);
    const auto *sig0_378 = buffer.data(sig0 + 378);
    const auto *sig0_380 = buffer.data(sig0 + 380);
    const auto *sig0_385 = buffer.data(sig0 + 385);
    const auto *sig0_387 = buffer.data(sig0 + 387);
    const auto *sig0_389 = buffer.data(sig0 + 389);
    const auto *sig0_393 = buffer.data(sig0 + 393);
    const auto *sig0_400 = buffer.data(sig0 + 400);
    const auto *sig0_402 = buffer.data(sig0 + 402);
    const auto *sig0_404 = buffer.data(sig0 + 404);
    const auto *sig0_405 = buffer.data(sig0 + 405);
    const auto *sig0_408 = buffer.data(sig0 + 408);
    const auto *sig0_410 = buffer.data(sig0 + 410);
    const auto *sig0_415 = buffer.data(sig0 + 415);
    const auto *sig0_417 = buffer.data(sig0 + 417);
    const auto *sig0_419 = buffer.data(sig0 + 419);

    const auto *sif_176 = buffer.data(sif + 176);
    const auto *sif_180 = buffer.data(sif + 180);
    const auto *sif_186 = buffer.data(sif + 186);
    const auto *sif_189 = buffer.data(sif + 189);
    const auto *sif_190 = buffer.data(sif + 190);
    const auto *sif_192 = buffer.data(sif + 192);
    const auto *sif_196 = buffer.data(sif + 196);
    const auto *sif_199 = buffer.data(sif + 199);
    const auto *sif_200 = buffer.data(sif + 200);
    const auto *sif_202 = buffer.data(sif + 202);
    const auto *sif_206 = buffer.data(sif + 206);
    const auto *sif_209 = buffer.data(sif + 209);
    const auto *sif_210 = buffer.data(sif + 210);
    const auto *sif_212 = buffer.data(sif + 212);
    const auto *sif_216 = buffer.data(sif + 216);
    const auto *sif_217 = buffer.data(sif + 217);
    const auto *sif_218 = buffer.data(sif + 218);
    const auto *sif_219 = buffer.data(sif + 219);
    const auto *sif_220 = buffer.data(sif + 220);
    const auto *sif_222 = buffer.data(sif + 222);
    const auto *sif_226 = buffer.data(sif + 226);
    const auto *sif_229 = buffer.data(sif + 229);
    const auto *sif_230 = buffer.data(sif + 230);
    const auto *sif_232 = buffer.data(sif + 232);
    const auto *sif_236 = buffer.data(sif + 236);
    const auto *sif_238 = buffer.data(sif + 238);
    const auto *sif_239 = buffer.data(sif + 239);
    const auto *sif_240 = buffer.data(sif + 240);
    const auto *sif_242 = buffer.data(sif + 242);
    const auto *sif_245 = buffer.data(sif + 245);
    const auto *sif_246 = buffer.data(sif + 246);
    const auto *sif_247 = buffer.data(sif + 247);
    const auto *sif_248 = buffer.data(sif + 248);
    const auto *sif_249 = buffer.data(sif + 249);
    const auto *sif_250 = buffer.data(sif + 250);
    const auto *sif_252 = buffer.data(sif + 252);
    const auto *sif_253 = buffer.data(sif + 253);
    const auto *sif_255 = buffer.data(sif + 255);
    const auto *sif_256 = buffer.data(sif + 256);
    const auto *sif_257 = buffer.data(sif + 257);
    const auto *sif_258 = buffer.data(sif + 258);
    const auto *sif_259 = buffer.data(sif + 259);
    const auto *sif_263 = buffer.data(sif + 263);
    const auto *sif_266 = buffer.data(sif + 266);
    const auto *sif_267 = buffer.data(sif + 267);
    const auto *sif_268 = buffer.data(sif + 268);
    const auto *sif_269 = buffer.data(sif + 269);
    const auto *sif_270 = buffer.data(sif + 270);
    const auto *sif_273 = buffer.data(sif + 273);
    const auto *sif_275 = buffer.data(sif + 275);
    const auto *sif_276 = buffer.data(sif + 276);
    const auto *sif_277 = buffer.data(sif + 277);
    const auto *sif_278 = buffer.data(sif + 278);
    const auto *sif_279 = buffer.data(sif + 279);

    const auto *sig1_300 = buffer.data(sig1 + 300);
    const auto *sig1_305 = buffer.data(sig1 + 305);
    const auto *sig1_315 = buffer.data(sig1 + 315);
    const auto *sig1_318 = buffer.data(sig1 + 318);
    const auto *sig1_325 = buffer.data(sig1 + 325);
    const auto *sig1_327 = buffer.data(sig1 + 327);
    const auto *sig1_365 = buffer.data(sig1 + 365);
    const auto *sig1_370 = buffer.data(sig1 + 370);
    const auto *sig1_372 = buffer.data(sig1 + 372);
    const auto *sig1_374 = buffer.data(sig1 + 374);
    const auto *sig1_375 = buffer.data(sig1 + 375);
    const auto *sig1_378 = buffer.data(sig1 + 378);
    const auto *sig1_380 = buffer.data(sig1 + 380);
    const auto *sig1_385 = buffer.data(sig1 + 385);
    const auto *sig1_387 = buffer.data(sig1 + 387);
    const auto *sig1_389 = buffer.data(sig1 + 389);
    const auto *sig1_393 = buffer.data(sig1 + 393);
    const auto *sig1_400 = buffer.data(sig1 + 400);
    const auto *sig1_402 = buffer.data(sig1 + 402);
    const auto *sig1_404 = buffer.data(sig1 + 404);
    const auto *sig1_405 = buffer.data(sig1 + 405);
    const auto *sig1_408 = buffer.data(sig1 + 408);
    const auto *sig1_410 = buffer.data(sig1 + 410);
    const auto *sig1_415 = buffer.data(sig1 + 415);
    const auto *sig1_417 = buffer.data(sig1 + 417);
    const auto *sig1_419 = buffer.data(sig1 + 419);

    const auto *skd0_168 = buffer.data(skd0 + 168);
    const auto *skd0_171 = buffer.data(skd0 + 171);
    const auto *skd0_173 = buffer.data(skd0 + 173);
    const auto *skd0_179 = buffer.data(skd0 + 179);
    const auto *skd0_180 = buffer.data(skd0 + 180);
    const auto *skd0_183 = buffer.data(skd0 + 183);
    const auto *skd0_185 = buffer.data(skd0 + 185);
    const auto *skd0_186 = buffer.data(skd0 + 186);
    const auto *skd0_189 = buffer.data(skd0 + 189);
    const auto *skd0_191 = buffer.data(skd0 + 191);
    const auto *skd0_192 = buffer.data(skd0 + 192);
    const auto *skd0_195 = buffer.data(skd0 + 195);
    const auto *skd0_197 = buffer.data(skd0 + 197);

    const auto *skd1_168 = buffer.data(skd1 + 168);
    const auto *skd1_171 = buffer.data(skd1 + 171);
    const auto *skd1_173 = buffer.data(skd1 + 173);
    const auto *skd1_179 = buffer.data(skd1 + 179);
    const auto *skd1_180 = buffer.data(skd1 + 180);
    const auto *skd1_183 = buffer.data(skd1 + 183);
    const auto *skd1_185 = buffer.data(skd1 + 185);
    const auto *skd1_186 = buffer.data(skd1 + 186);
    const auto *skd1_189 = buffer.data(skd1 + 189);
    const auto *skd1_191 = buffer.data(skd1 + 191);
    const auto *skd1_192 = buffer.data(skd1 + 192);
    const auto *skd1_195 = buffer.data(skd1 + 195);
    const auto *skd1_197 = buffer.data(skd1 + 197);

    const auto *skf_246 = buffer.data(skf + 246);
    const auto *skf_247 = buffer.data(skf + 247);
    const auto *skf_248 = buffer.data(skf + 248);
    const auto *skf_249 = buffer.data(skf + 249);
    const auto *skf_250 = buffer.data(skf + 250);
    const auto *skf_252 = buffer.data(skf + 252);
    const auto *skf_256 = buffer.data(skf + 256);
    const auto *skf_257 = buffer.data(skf + 257);
    const auto *skf_258 = buffer.data(skf + 258);
    const auto *skf_259 = buffer.data(skf + 259);
    const auto *skf_260 = buffer.data(skf + 260);
    const auto *skf_262 = buffer.data(skf + 262);
    const auto *skf_266 = buffer.data(skf + 266);
    const auto *skf_267 = buffer.data(skf + 267);
    const auto *skf_268 = buffer.data(skf + 268);
    const auto *skf_269 = buffer.data(skf + 269);
    const auto *skf_270 = buffer.data(skf + 270);
    const auto *skf_272 = buffer.data(skf + 272);
    const auto *skf_276 = buffer.data(skf + 276);
    const auto *skf_277 = buffer.data(skf + 277);
    const auto *skf_278 = buffer.data(skf + 278);
    const auto *skf_279 = buffer.data(skf + 279);
    const auto *skf_280 = buffer.data(skf + 280);
    const auto *skf_282 = buffer.data(skf + 282);
    const auto *skf_283 = buffer.data(skf + 283);
    const auto *skf_285 = buffer.data(skf + 285);
    const auto *skf_286 = buffer.data(skf + 286);
    const auto *skf_287 = buffer.data(skf + 287);
    const auto *skf_288 = buffer.data(skf + 288);
    const auto *skf_289 = buffer.data(skf + 289);
    const auto *skf_290 = buffer.data(skf + 290);
    const auto *skf_292 = buffer.data(skf + 292);
    const auto *skf_295 = buffer.data(skf + 295);
    const auto *skf_296 = buffer.data(skf + 296);
    const auto *skf_297 = buffer.data(skf + 297);
    const auto *skf_298 = buffer.data(skf + 298);
    const auto *skf_299 = buffer.data(skf + 299);
    const auto *skf_300 = buffer.data(skf + 300);
    const auto *skf_302 = buffer.data(skf + 302);
    const auto *skf_303 = buffer.data(skf + 303);
    const auto *skf_305 = buffer.data(skf + 305);
    const auto *skf_306 = buffer.data(skf + 306);
    const auto *skf_307 = buffer.data(skf + 307);
    const auto *skf_308 = buffer.data(skf + 308);
    const auto *skf_309 = buffer.data(skf + 309);
    const auto *skf_310 = buffer.data(skf + 310);
    const auto *skf_312 = buffer.data(skf + 312);
    const auto *skf_313 = buffer.data(skf + 313);
    const auto *skf_315 = buffer.data(skf + 315);
    const auto *skf_316 = buffer.data(skf + 316);
    const auto *skf_317 = buffer.data(skf + 317);
    const auto *skf_318 = buffer.data(skf + 318);
    const auto *skf_319 = buffer.data(skf + 319);
    const auto *skf_320 = buffer.data(skf + 320);
    const auto *skf_322 = buffer.data(skf + 322);
    const auto *skf_323 = buffer.data(skf + 323);
    const auto *skf_325 = buffer.data(skf + 325);
    const auto *skf_326 = buffer.data(skf + 326);
    const auto *skf_327 = buffer.data(skf + 327);
    const auto *skf_328 = buffer.data(skf + 328);
    const auto *skf_329 = buffer.data(skf + 329);

#pragma omp simd aligned(t_365, t_366, t_367, t_368, pb_x, pc_x, sig0_365, sif_245, sif_246, \
                         sif_247, sif_248, sig1_365, skf_246, skf_247, \
                         skf_248 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_365[k] = pb_x[k] * sig0_365[k]
                   + f_8 * sif_245[k]
                   - f_6 * pc_x[k] * sig1_365[k];

        t_366[k] = f_7 * sif_246[k]
                   + f_3 * pc_x[k] * skf_246[k];

        t_367[k] = f_7 * sif_247[k]
                   + f_3 * pc_x[k] * skf_247[k];

        t_368[k] = f_7 * sif_248[k]
                   + f_3 * pc_x[k] * skf_248[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, t_372, pb_x, pc_x, pc_z, sig0_370, sig0_372, \
                         sif_176, sif_249, sig1_370, sig1_372, skf_246, \
                         skf_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = f_7 * sif_249[k]
                   + f_3 * pc_x[k] * skf_249[k];

        t_370[k] = pb_x[k] * sig0_370[k]
                   - f_6 * pc_x[k] * sig1_370[k];

        t_371[k] = f_12 * sif_176[k]
                   + f_3 * pc_z[k] * skf_246[k];

        t_372[k] = pb_x[k] * sig0_372[k]
                   - f_6 * pc_x[k] * sig1_372[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, t_376, pb_x, pc_x, pc_y, sig0_374, sig0_375, \
                         sif_189, sif_190, sif_250, sig1_374, sig1_375, skf_249, \
                         skf_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_12 * sif_189[k]
                   + f_3 * pc_y[k] * skf_249[k];

        t_374[k] = pb_x[k] * sig0_374[k]
                   - f_6 * pc_x[k] * sig1_374[k];

        t_375[k] = pb_x[k] * sig0_375[k]
                   + f_11 * sif_250[k]
                   - f_6 * pc_x[k] * sig1_375[k];

        t_376[k] = f_8 * sif_190[k]
                   + f_3 * pc_y[k] * skf_250[k];
    }

#pragma omp simd aligned(t_377, t_378, t_379, pb_x, pc_x, pc_y, pc_z, sig0_378, sif_180, \
                         sif_192, sif_253, sig1_378, skf_250, skf_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = f_11 * sif_180[k]
                   + f_3 * pc_z[k] * skf_250[k];

        t_378[k] = pb_x[k] * sig0_378[k]
                   + f_8 * sif_253[k]
                   - f_6 * pc_x[k] * sig1_378[k];

        t_379[k] = f_8 * sif_192[k]
                   + f_3 * pc_y[k] * skf_252[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, pb_x, pc_x, sig0_380, sif_255, sif_256, \
                         sif_257, sif_258, sig1_380, skf_256, skf_257, \
                         skf_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = pb_x[k] * sig0_380[k]
                   + f_8 * sif_255[k]
                   - f_6 * pc_x[k] * sig1_380[k];

        t_381[k] = f_7 * sif_256[k]
                   + f_3 * pc_x[k] * skf_256[k];

        t_382[k] = f_7 * sif_257[k]
                   + f_3 * pc_x[k] * skf_257[k];

        t_383[k] = f_7 * sif_258[k]
                   + f_3 * pc_x[k] * skf_258[k];
    }

#pragma omp simd aligned(t_384, t_385, t_386, t_387, pb_x, pc_x, pc_z, sig0_385, sig0_387, \
                         sif_186, sif_259, sig1_385, sig1_387, skf_256, \
                         skf_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_384[k] = f_7 * sif_259[k]
                   + f_3 * pc_x[k] * skf_259[k];

        t_385[k] = pb_x[k] * sig0_385[k]
                   - f_6 * pc_x[k] * sig1_385[k];

        t_386[k] = f_11 * sif_186[k]
                   + f_3 * pc_z[k] * skf_256[k];

        t_387[k] = pb_x[k] * sig0_387[k]
                   - f_6 * pc_x[k] * sig1_387[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, t_391, pb_x, pb_y, pc_x, pc_y, sig0_300, \
                         sig0_389, sif_199, sif_200, sig1_300, sig1_389, skf_259, \
                         skf_260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_8 * sif_199[k]
                   + f_3 * pc_y[k] * skf_259[k];

        t_389[k] = pb_x[k] * sig0_389[k]
                   - f_6 * pc_x[k] * sig1_389[k];

        t_390[k] = pb_y[k] * sig0_300[k]
                   - f_6 * pc_y[k] * sig1_300[k];

        t_391[k] = f_7 * sif_200[k]
                   + f_3 * pc_y[k] * skf_260[k];
    }

#pragma omp simd aligned(t_392, t_393, t_394, pb_x, pc_x, pc_y, pc_z, sig0_393, sif_190, \
                         sif_202, sif_263, sig1_393, skf_260, skf_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_392[k] = f_10 * sif_190[k]
                   + f_3 * pc_z[k] * skf_260[k];

        t_393[k] = pb_x[k] * sig0_393[k]
                   + f_8 * sif_263[k]
                   - f_6 * pc_x[k] * sig1_393[k];

        t_394[k] = f_7 * sif_202[k]
                   + f_3 * pc_y[k] * skf_262[k];
    }

#pragma omp simd aligned(t_395, t_396, t_397, t_398, pb_y, pc_x, pc_y, sig0_305, sif_266, \
                         sif_267, sif_268, sig1_305, skf_266, skf_267, \
                         skf_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = pb_y[k] * sig0_305[k]
                   - f_6 * pc_y[k] * sig1_305[k];

        t_396[k] = f_7 * sif_266[k]
                   + f_3 * pc_x[k] * skf_266[k];

        t_397[k] = f_7 * sif_267[k]
                   + f_3 * pc_x[k] * skf_267[k];

        t_398[k] = f_7 * sif_268[k]
                   + f_3 * pc_x[k] * skf_268[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, t_402, pb_x, pc_x, pc_z, sig0_400, sig0_402, \
                         sif_196, sif_269, sig1_400, sig1_402, skf_266, \
                         skf_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_399[k] = f_7 * sif_269[k]
                   + f_3 * pc_x[k] * skf_269[k];

        t_400[k] = pb_x[k] * sig0_400[k]
                   - f_6 * pc_x[k] * sig1_400[k];

        t_401[k] = f_10 * sif_196[k]
                   + f_3 * pc_z[k] * skf_266[k];

        t_402[k] = pb_x[k] * sig0_402[k]
                   - f_6 * pc_x[k] * sig1_402[k];
    }

#pragma omp simd aligned(t_403, t_404, t_405, t_406, pb_x, pc_x, pc_y, sig0_404, sig0_405, \
                         sif_209, sif_270, sig1_404, sig1_405, skf_269, \
                         skf_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_403[k] = f_7 * sif_209[k]
                   + f_3 * pc_y[k] * skf_269[k];

        t_404[k] = pb_x[k] * sig0_404[k]
                   - f_6 * pc_x[k] * sig1_404[k];

        t_405[k] = pb_x[k] * sig0_405[k]
                   + f_11 * sif_270[k]
                   - f_6 * pc_x[k] * sig1_405[k];

        t_406[k] = f_3 * pc_y[k] * skf_270[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, pb_x, pc_x, pc_y, pc_z, sig0_408, sif_200, \
                         sif_273, sig1_408, skf_270, skf_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = f_9 * sif_200[k]
                   + f_3 * pc_z[k] * skf_270[k];

        t_408[k] = pb_x[k] * sig0_408[k]
                   + f_8 * sif_273[k]
                   - f_6 * pc_x[k] * sig1_408[k];

        t_409[k] = f_3 * pc_y[k] * skf_272[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, pb_x, pc_x, sig0_410, sif_275, sif_276, \
                         sif_277, sif_278, sig1_410, skf_276, skf_277, \
                         skf_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = pb_x[k] * sig0_410[k]
                   + f_8 * sif_275[k]
                   - f_6 * pc_x[k] * sig1_410[k];

        t_411[k] = f_7 * sif_276[k]
                   + f_3 * pc_x[k] * skf_276[k];

        t_412[k] = f_7 * sif_277[k]
                   + f_3 * pc_x[k] * skf_277[k];

        t_413[k] = f_7 * sif_278[k]
                   + f_3 * pc_x[k] * skf_278[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, pb_x, pc_x, pc_z, sig0_415, sig0_417, \
                         sif_206, sif_279, sig1_415, sig1_417, skf_276, \
                         skf_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_7 * sif_279[k]
                   + f_3 * pc_x[k] * skf_279[k];

        t_415[k] = pb_x[k] * sig0_415[k]
                   - f_6 * pc_x[k] * sig1_415[k];

        t_416[k] = f_9 * sif_206[k]
                   + f_3 * pc_z[k] * skf_276[k];

        t_417[k] = pb_x[k] * sig0_417[k]
                   - f_6 * pc_x[k] * sig1_417[k];
    }

#pragma omp simd aligned(t_418, t_419, t_420, t_421, t_422, pb_x, pc_x, pc_y, pc_z, sig0_419, \
                         sif_210, sig1_419, skd0_168, skd1_168, skf_279, \
                         skf_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_418[k] = f_3 * pc_y[k] * skf_279[k];

        t_419[k] = pb_x[k] * sig0_419[k]
                   - f_6 * pc_x[k] * sig1_419[k];

        t_420[k] = f_1 * skd0_168[k]
                   - f_2 * skd1_168[k]
                   + f_3 * pc_x[k] * skf_280[k];

        t_421[k] = f_0 * sif_210[k]
                   + f_3 * pc_y[k] * skf_280[k];

        t_422[k] = f_3 * pc_z[k] * skf_280[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, t_426, pc_x, pc_y, sif_212, skd0_171, skd0_173, \
                         skd1_171, skd1_173, skf_282, skf_283, skf_285, \
                         skf_286 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = f_4 * skd0_171[k]
                   - f_5 * skd1_171[k]
                   + f_3 * pc_x[k] * skf_283[k];

        t_424[k] = f_0 * sif_212[k]
                   + f_3 * pc_y[k] * skf_282[k];

        t_425[k] = f_4 * skd0_173[k]
                   - f_5 * skd1_173[k]
                   + f_3 * pc_x[k] * skf_285[k];

        t_426[k] = f_3 * pc_x[k] * skf_286[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, t_430, t_431, pc_x, pc_y, pc_z, sif_216, \
                         skd0_171, skd1_171, skf_286, skf_287, skf_288, \
                         skf_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = f_3 * pc_x[k] * skf_287[k];

        t_428[k] = f_3 * pc_x[k] * skf_288[k];

        t_429[k] = f_3 * pc_x[k] * skf_289[k];

        t_430[k] = f_0 * sif_216[k]
                   + f_1 * skd0_171[k]
                   - f_2 * skd1_171[k]
                   + f_3 * pc_y[k] * skf_286[k];

        t_431[k] = f_3 * pc_z[k] * skf_286[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, pb_z, pc_y, pc_z, sig0_315, sif_218, \
                         sif_219, sig1_315, skd0_173, skd1_173, skf_288, \
                         skf_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = f_0 * sif_218[k]
                   + f_4 * skd0_173[k]
                   - f_5 * skd1_173[k]
                   + f_3 * pc_y[k] * skf_288[k];

        t_433[k] = f_0 * sif_219[k]
                   + f_3 * pc_y[k] * skf_289[k];

        t_434[k] = f_1 * skd0_173[k]
                   - f_2 * skd1_173[k]
                   + f_3 * pc_z[k] * skf_289[k];

        t_435[k] = pb_z[k] * sig0_315[k]
                   - f_6 * pc_z[k] * sig1_315[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, pb_z, pc_y, pc_z, sig0_318, sif_210, \
                         sif_220, sif_222, sig1_318, skf_290, skf_292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_9 * sif_220[k]
                   + f_3 * pc_y[k] * skf_290[k];

        t_437[k] = f_7 * sif_210[k]
                   + f_3 * pc_z[k] * skf_290[k];

        t_438[k] = pb_z[k] * sig0_318[k]
                   - f_6 * pc_z[k] * sig1_318[k];

        t_439[k] = f_9 * sif_222[k]
                   + f_3 * pc_y[k] * skf_292[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, pc_x, skd0_179, skd1_179, skf_295, \
                         skf_296, skf_297, skf_298, skf_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_4 * skd0_179[k]
                   - f_5 * skd1_179[k]
                   + f_3 * pc_x[k] * skf_295[k];

        t_441[k] = f_3 * pc_x[k] * skf_296[k];

        t_442[k] = f_3 * pc_x[k] * skf_297[k];

        t_443[k] = f_3 * pc_x[k] * skf_298[k];

        t_444[k] = f_3 * pc_x[k] * skf_299[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, pb_z, pc_y, pc_z, sig0_325, sig0_327, \
                         sif_216, sif_217, sif_229, sig1_325, sig1_327, skf_296, \
                         skf_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = pb_z[k] * sig0_325[k]
                   - f_6 * pc_z[k] * sig1_325[k];

        t_446[k] = f_7 * sif_216[k]
                   + f_3 * pc_z[k] * skf_296[k];

        t_447[k] = pb_z[k] * sig0_327[k]
                   + f_8 * sif_217[k]
                   - f_6 * pc_z[k] * sig1_327[k];

        t_448[k] = f_9 * sif_229[k]
                   + f_3 * pc_y[k] * skf_299[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, pc_x, pc_y, pc_z, sif_219, sif_220, \
                         sif_230, skd0_179, skd0_180, skd1_179, skd1_180, skf_299, \
                         skf_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_7 * sif_219[k]
                   + f_1 * skd0_179[k]
                   - f_2 * skd1_179[k]
                   + f_3 * pc_z[k] * skf_299[k];

        t_450[k] = f_1 * skd0_180[k]
                   - f_2 * skd1_180[k]
                   + f_3 * pc_x[k] * skf_300[k];

        t_451[k] = f_10 * sif_230[k]
                   + f_3 * pc_y[k] * skf_300[k];

        t_452[k] = f_8 * sif_220[k]
                   + f_3 * pc_z[k] * skf_300[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, t_456, pc_x, pc_y, sif_232, skd0_183, skd0_185, \
                         skd1_183, skd1_185, skf_302, skf_303, skf_305, \
                         skf_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = f_4 * skd0_183[k]
                   - f_5 * skd1_183[k]
                   + f_3 * pc_x[k] * skf_303[k];

        t_454[k] = f_10 * sif_232[k]
                   + f_3 * pc_y[k] * skf_302[k];

        t_455[k] = f_4 * skd0_185[k]
                   - f_5 * skd1_185[k]
                   + f_3 * pc_x[k] * skf_305[k];

        t_456[k] = f_3 * pc_x[k] * skf_306[k];
    }

#pragma omp simd aligned(t_457, t_458, t_459, t_460, t_461, pc_x, pc_y, pc_z, sif_226, \
                         sif_236, skd0_183, skd1_183, skf_306, skf_307, skf_308, \
                         skf_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_457[k] = f_3 * pc_x[k] * skf_307[k];

        t_458[k] = f_3 * pc_x[k] * skf_308[k];

        t_459[k] = f_3 * pc_x[k] * skf_309[k];

        t_460[k] = f_10 * sif_236[k]
                   + f_1 * skd0_183[k]
                   - f_2 * skd1_183[k]
                   + f_3 * pc_y[k] * skf_306[k];

        t_461[k] = f_8 * sif_226[k]
                   + f_3 * pc_z[k] * skf_306[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, pc_y, pc_z, sif_229, sif_238, sif_239, skd0_185, \
                         skd1_185, skf_308, skf_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = f_10 * sif_238[k]
                   + f_4 * skd0_185[k]
                   - f_5 * skd1_185[k]
                   + f_3 * pc_y[k] * skf_308[k];

        t_463[k] = f_10 * sif_239[k]
                   + f_3 * pc_y[k] * skf_309[k];

        t_464[k] = f_8 * sif_229[k]
                   + f_1 * skd0_185[k]
                   - f_2 * skd1_185[k]
                   + f_3 * pc_z[k] * skf_309[k];
    }

#pragma omp simd aligned(t_465, t_466, t_467, t_468, pc_x, pc_y, pc_z, sif_230, sif_240, \
                         skd0_186, skd0_189, skd1_186, skd1_189, skf_310, \
                         skf_313 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_465[k] = f_1 * skd0_186[k]
                   - f_2 * skd1_186[k]
                   + f_3 * pc_x[k] * skf_310[k];

        t_466[k] = f_11 * sif_240[k]
                   + f_3 * pc_y[k] * skf_310[k];

        t_467[k] = f_12 * sif_230[k]
                   + f_3 * pc_z[k] * skf_310[k];

        t_468[k] = f_4 * skd0_189[k]
                   - f_5 * skd1_189[k]
                   + f_3 * pc_x[k] * skf_313[k];
    }

#pragma omp simd aligned(t_469, t_470, t_471, t_472, t_473, pc_x, pc_y, sif_242, skd0_191, \
                         skd1_191, skf_312, skf_315, skf_316, skf_317, \
                         skf_318 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_469[k] = f_11 * sif_242[k]
                   + f_3 * pc_y[k] * skf_312[k];

        t_470[k] = f_4 * skd0_191[k]
                   - f_5 * skd1_191[k]
                   + f_3 * pc_x[k] * skf_315[k];

        t_471[k] = f_3 * pc_x[k] * skf_316[k];

        t_472[k] = f_3 * pc_x[k] * skf_317[k];

        t_473[k] = f_3 * pc_x[k] * skf_318[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, pc_x, pc_y, pc_z, sif_236, sif_246, skd0_189, \
                         skd1_189, skf_316, skf_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_3 * pc_x[k] * skf_319[k];

        t_475[k] = f_11 * sif_246[k]
                   + f_1 * skd0_189[k]
                   - f_2 * skd1_189[k]
                   + f_3 * pc_y[k] * skf_316[k];

        t_476[k] = f_12 * sif_236[k]
                   + f_3 * pc_z[k] * skf_316[k];
    }

#pragma omp simd aligned(t_477, t_478, t_479, pc_y, pc_z, sif_239, sif_248, sif_249, skd0_191, \
                         skd1_191, skf_318, skf_319 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_477[k] = f_11 * sif_248[k]
                   + f_4 * skd0_191[k]
                   - f_5 * skd1_191[k]
                   + f_3 * pc_y[k] * skf_318[k];

        t_478[k] = f_11 * sif_249[k]
                   + f_3 * pc_y[k] * skf_319[k];

        t_479[k] = f_12 * sif_239[k]
                   + f_1 * skd0_191[k]
                   - f_2 * skd1_191[k]
                   + f_3 * pc_z[k] * skf_319[k];
    }

#pragma omp simd aligned(t_480, t_481, t_482, t_483, pc_x, pc_y, pc_z, sif_240, sif_250, \
                         skd0_192, skd0_195, skd1_192, skd1_195, skf_320, \
                         skf_323 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_480[k] = f_1 * skd0_192[k]
                   - f_2 * skd1_192[k]
                   + f_3 * pc_x[k] * skf_320[k];

        t_481[k] = f_12 * sif_250[k]
                   + f_3 * pc_y[k] * skf_320[k];

        t_482[k] = f_11 * sif_240[k]
                   + f_3 * pc_z[k] * skf_320[k];

        t_483[k] = f_4 * skd0_195[k]
                   - f_5 * skd1_195[k]
                   + f_3 * pc_x[k] * skf_323[k];
    }

#pragma omp simd aligned(t_484, t_485, t_486, t_487, t_488, pc_x, pc_y, sif_252, skd0_197, \
                         skd1_197, skf_322, skf_325, skf_326, skf_327, \
                         skf_328 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_484[k] = f_12 * sif_252[k]
                   + f_3 * pc_y[k] * skf_322[k];

        t_485[k] = f_4 * skd0_197[k]
                   - f_5 * skd1_197[k]
                   + f_3 * pc_x[k] * skf_325[k];

        t_486[k] = f_3 * pc_x[k] * skf_326[k];

        t_487[k] = f_3 * pc_x[k] * skf_327[k];

        t_488[k] = f_3 * pc_x[k] * skf_328[k];
    }

#pragma omp simd aligned(t_489, t_490, t_491, pc_x, pc_y, pc_z, sif_246, sif_256, skd0_195, \
                         skd1_195, skf_326, skf_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_489[k] = f_3 * pc_x[k] * skf_329[k];

        t_490[k] = f_12 * sif_256[k]
                   + f_1 * skd0_195[k]
                   - f_2 * skd1_195[k]
                   + f_3 * pc_y[k] * skf_326[k];

        t_491[k] = f_11 * sif_246[k]
                   + f_3 * pc_z[k] * skf_326[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, pc_y, pc_z, sif_249, sif_258, sif_259, skd0_197, \
                         skd1_197, skf_328, skf_329 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = f_12 * sif_258[k]
                   + f_4 * skd0_197[k]
                   - f_5 * skd1_197[k]
                   + f_3 * pc_y[k] * skf_328[k];

        t_493[k] = f_12 * sif_259[k]
                   + f_3 * pc_y[k] * skf_329[k];

        t_494[k] = f_11 * sif_249[k]
                   + f_1 * skd0_197[k]
                   - f_2 * skd1_197[k]
                   + f_3 * pc_z[k] * skf_329[k];
    }
}

static auto
compute_prim_skg_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pb,
                                                          const size_t pc, const size_t sig0,
                                                          const size_t sif, const size_t sig1,
                                                          const size_t skd0, const size_t skd1,
                                                          const size_t skf, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 3.0 / q;
    const auto f_10 = 2.5 / q;
    const auto f_11 = 2.0 / q;

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

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sig0_405 = buffer.data(sig0 + 405);
    const auto *sig0_410 = buffer.data(sig0 + 410);
    const auto *sig0_415 = buffer.data(sig0 + 415);
    const auto *sig0_417 = buffer.data(sig0 + 417);
    const auto *sig0_419 = buffer.data(sig0 + 419);

    const auto *sif_250 = buffer.data(sif + 250);
    const auto *sif_256 = buffer.data(sif + 256);
    const auto *sif_259 = buffer.data(sif + 259);
    const auto *sif_260 = buffer.data(sif + 260);
    const auto *sif_262 = buffer.data(sif + 262);
    const auto *sif_266 = buffer.data(sif + 266);
    const auto *sif_268 = buffer.data(sif + 268);
    const auto *sif_269 = buffer.data(sif + 269);
    const auto *sif_270 = buffer.data(sif + 270);
    const auto *sif_272 = buffer.data(sif + 272);
    const auto *sif_276 = buffer.data(sif + 276);
    const auto *sif_278 = buffer.data(sif + 278);
    const auto *sif_279 = buffer.data(sif + 279);

    const auto *sig1_405 = buffer.data(sig1 + 405);
    const auto *sig1_410 = buffer.data(sig1 + 410);
    const auto *sig1_415 = buffer.data(sig1 + 415);
    const auto *sig1_417 = buffer.data(sig1 + 417);
    const auto *sig1_419 = buffer.data(sig1 + 419);

    const auto *skd0_198 = buffer.data(skd0 + 198);
    const auto *skd0_201 = buffer.data(skd0 + 201);
    const auto *skd0_203 = buffer.data(skd0 + 203);
    const auto *skd0_207 = buffer.data(skd0 + 207);
    const auto *skd0_210 = buffer.data(skd0 + 210);
    const auto *skd0_213 = buffer.data(skd0 + 213);
    const auto *skd0_215 = buffer.data(skd0 + 215);

    const auto *skd1_198 = buffer.data(skd1 + 198);
    const auto *skd1_201 = buffer.data(skd1 + 201);
    const auto *skd1_203 = buffer.data(skd1 + 203);
    const auto *skd1_207 = buffer.data(skd1 + 207);
    const auto *skd1_210 = buffer.data(skd1 + 210);
    const auto *skd1_213 = buffer.data(skd1 + 213);
    const auto *skd1_215 = buffer.data(skd1 + 215);

    const auto *skf_330 = buffer.data(skf + 330);
    const auto *skf_332 = buffer.data(skf + 332);
    const auto *skf_333 = buffer.data(skf + 333);
    const auto *skf_335 = buffer.data(skf + 335);
    const auto *skf_336 = buffer.data(skf + 336);
    const auto *skf_337 = buffer.data(skf + 337);
    const auto *skf_338 = buffer.data(skf + 338);
    const auto *skf_339 = buffer.data(skf + 339);
    const auto *skf_340 = buffer.data(skf + 340);
    const auto *skf_342 = buffer.data(skf + 342);
    const auto *skf_343 = buffer.data(skf + 343);
    const auto *skf_346 = buffer.data(skf + 346);
    const auto *skf_347 = buffer.data(skf + 347);
    const auto *skf_348 = buffer.data(skf + 348);
    const auto *skf_349 = buffer.data(skf + 349);
    const auto *skf_350 = buffer.data(skf + 350);
    const auto *skf_352 = buffer.data(skf + 352);
    const auto *skf_353 = buffer.data(skf + 353);
    const auto *skf_355 = buffer.data(skf + 355);
    const auto *skf_356 = buffer.data(skf + 356);
    const auto *skf_357 = buffer.data(skf + 357);
    const auto *skf_358 = buffer.data(skf + 358);
    const auto *skf_359 = buffer.data(skf + 359);

#pragma omp simd aligned(t_495, t_496, t_497, t_498, pc_x, pc_y, pc_z, sif_250, sif_260, \
                         skd0_198, skd0_201, skd1_198, skd1_201, skf_330, \
                         skf_333 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_495[k] = f_1 * skd0_198[k]
                   - f_2 * skd1_198[k]
                   + f_3 * pc_x[k] * skf_330[k];

        t_496[k] = f_8 * sif_260[k]
                   + f_3 * pc_y[k] * skf_330[k];

        t_497[k] = f_10 * sif_250[k]
                   + f_3 * pc_z[k] * skf_330[k];

        t_498[k] = f_4 * skd0_201[k]
                   - f_5 * skd1_201[k]
                   + f_3 * pc_x[k] * skf_333[k];
    }

#pragma omp simd aligned(t_499, t_500, t_501, t_502, t_503, pc_x, pc_y, sif_262, skd0_203, \
                         skd1_203, skf_332, skf_335, skf_336, skf_337, \
                         skf_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_499[k] = f_8 * sif_262[k]
                   + f_3 * pc_y[k] * skf_332[k];

        t_500[k] = f_4 * skd0_203[k]
                   - f_5 * skd1_203[k]
                   + f_3 * pc_x[k] * skf_335[k];

        t_501[k] = f_3 * pc_x[k] * skf_336[k];

        t_502[k] = f_3 * pc_x[k] * skf_337[k];

        t_503[k] = f_3 * pc_x[k] * skf_338[k];
    }

#pragma omp simd aligned(t_504, t_505, t_506, pc_x, pc_y, pc_z, sif_256, sif_266, skd0_201, \
                         skd1_201, skf_336, skf_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_504[k] = f_3 * pc_x[k] * skf_339[k];

        t_505[k] = f_8 * sif_266[k]
                   + f_1 * skd0_201[k]
                   - f_2 * skd1_201[k]
                   + f_3 * pc_y[k] * skf_336[k];

        t_506[k] = f_10 * sif_256[k]
                   + f_3 * pc_z[k] * skf_336[k];
    }

#pragma omp simd aligned(t_507, t_508, t_509, t_510, pb_y, pc_y, pc_z, sig0_405, sif_259, \
                         sif_268, sif_269, sig1_405, skd0_203, skd1_203, skf_338, \
                         skf_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_507[k] = f_8 * sif_268[k]
                   + f_4 * skd0_203[k]
                   - f_5 * skd1_203[k]
                   + f_3 * pc_y[k] * skf_338[k];

        t_508[k] = f_8 * sif_269[k]
                   + f_3 * pc_y[k] * skf_339[k];

        t_509[k] = f_10 * sif_259[k]
                   + f_1 * skd0_203[k]
                   - f_2 * skd1_203[k]
                   + f_3 * pc_z[k] * skf_339[k];

        t_510[k] = pb_y[k] * sig0_405[k]
                   - f_6 * pc_y[k] * sig1_405[k];
    }

#pragma omp simd aligned(t_511, t_512, t_513, t_514, pc_x, pc_y, pc_z, sif_260, sif_270, \
                         sif_272, skd0_207, skd1_207, skf_340, skf_342, \
                         skf_343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_511[k] = f_7 * sif_270[k]
                   + f_3 * pc_y[k] * skf_340[k];

        t_512[k] = f_9 * sif_260[k]
                   + f_3 * pc_z[k] * skf_340[k];

        t_513[k] = f_4 * skd0_207[k]
                   - f_5 * skd1_207[k]
                   + f_3 * pc_x[k] * skf_343[k];

        t_514[k] = f_7 * sif_272[k]
                   + f_3 * pc_y[k] * skf_342[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, t_518, t_519, pb_y, pc_x, pc_y, sig0_410, \
                         sig1_410, skf_346, skf_347, skf_348, skf_349 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = pb_y[k] * sig0_410[k]
                   - f_6 * pc_y[k] * sig1_410[k];

        t_516[k] = f_3 * pc_x[k] * skf_346[k];

        t_517[k] = f_3 * pc_x[k] * skf_347[k];

        t_518[k] = f_3 * pc_x[k] * skf_348[k];

        t_519[k] = f_3 * pc_x[k] * skf_349[k];
    }

#pragma omp simd aligned(t_520, t_521, t_522, pb_y, pc_y, pc_z, sig0_415, sig0_417, sif_266, \
                         sif_276, sif_278, sig1_415, sig1_417, \
                         skf_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_520[k] = pb_y[k] * sig0_415[k]
                   + f_11 * sif_276[k]
                   - f_6 * pc_y[k] * sig1_415[k];

        t_521[k] = f_9 * sif_266[k]
                   + f_3 * pc_z[k] * skf_346[k];

        t_522[k] = pb_y[k] * sig0_417[k]
                   + f_8 * sif_278[k]
                   - f_6 * pc_y[k] * sig1_417[k];
    }

#pragma omp simd aligned(t_523, t_524, t_525, t_526, pb_y, pc_x, pc_y, sig0_419, sif_279, \
                         sig1_419, skd0_210, skd1_210, skf_349, \
                         skf_350 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_523[k] = f_7 * sif_279[k]
                   + f_3 * pc_y[k] * skf_349[k];

        t_524[k] = pb_y[k] * sig0_419[k]
                   - f_6 * pc_y[k] * sig1_419[k];

        t_525[k] = f_1 * skd0_210[k]
                   - f_2 * skd1_210[k]
                   + f_3 * pc_x[k] * skf_350[k];

        t_526[k] = f_3 * pc_y[k] * skf_350[k];
    }

#pragma omp simd aligned(t_527, t_528, t_529, t_530, pc_x, pc_y, pc_z, sif_270, skd0_213, \
                         skd0_215, skd1_213, skd1_215, skf_350, skf_352, skf_353, \
                         skf_355 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_527[k] = f_0 * sif_270[k]
                   + f_3 * pc_z[k] * skf_350[k];

        t_528[k] = f_4 * skd0_213[k]
                   - f_5 * skd1_213[k]
                   + f_3 * pc_x[k] * skf_353[k];

        t_529[k] = f_3 * pc_y[k] * skf_352[k];

        t_530[k] = f_4 * skd0_215[k]
                   - f_5 * skd1_215[k]
                   + f_3 * pc_x[k] * skf_355[k];
    }

#pragma omp simd aligned(t_531, t_532, t_533, t_534, t_535, t_536, pc_x, pc_y, pc_z, sif_276, \
                         skd0_213, skd1_213, skf_356, skf_357, skf_358, \
                         skf_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_531[k] = f_3 * pc_x[k] * skf_356[k];

        t_532[k] = f_3 * pc_x[k] * skf_357[k];

        t_533[k] = f_3 * pc_x[k] * skf_358[k];

        t_534[k] = f_3 * pc_x[k] * skf_359[k];

        t_535[k] = f_1 * skd0_213[k]
                   - f_2 * skd1_213[k]
                   + f_3 * pc_y[k] * skf_356[k];

        t_536[k] = f_0 * sif_276[k]
                   + f_3 * pc_z[k] * skf_356[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, pc_y, pc_z, sif_279, skd0_215, skd1_215, \
                         skf_358, skf_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = f_4 * skd0_215[k]
                   - f_5 * skd1_215[k]
                   + f_3 * pc_y[k] * skf_358[k];

        t_538[k] = f_3 * pc_y[k] * skf_359[k];

        t_539[k] = f_0 * sif_279[k]
                   + f_1 * skd0_215[k]
                   - f_2 * skd1_215[k]
                   + f_3 * pc_z[k] * skf_359[k];
    }
}

auto
compute_prim_skg_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sig0, const size_t sif,
                                                   const size_t sig1, const size_t skd0,
                                                   const size_t skd1, const size_t skf,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_skg_three_center_electron_repulsion_0_piece0(buffer, target, pb, pc, sig0, sif,
                                                              sig1, skd0, skd1, skf, ncols,
                                                              gamma, p, q);

    compute_prim_skg_three_center_electron_repulsion_0_piece1(buffer, target, pb, pc, sig0, sif,
                                                              sig1, skd0, skd1, skf, ncols,
                                                              gamma, p, q);

    compute_prim_skg_three_center_electron_repulsion_0_piece2(buffer, target, pb, pc, sig0, sif,
                                                              sig1, skd0, skd1, skf, ncols,
                                                              gamma, p, q);

    compute_prim_skg_three_center_electron_repulsion_0_piece3(buffer, target, pb, pc, sig0, sif,
                                                              sig1, skd0, skd1, skf, ncols,
                                                              gamma, p, q);

    compute_prim_skg_three_center_electron_repulsion_0_piece4(buffer, target, pb, pc, sig0, sif,
                                                              sig1, skd0, skd1, skf, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
