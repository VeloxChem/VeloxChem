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


#include "SimdThreeCenterElectronRepulsionVrrRecSDH.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_sdh_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sph0, const size_t spg,
                                                   const size_t sph1, const size_t sdf0,
                                                   const size_t sdf1, const size_t sdg,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
    const auto f_1 = 2.0 / gamma;
    const auto f_2 = 2.0 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 1.0 / gamma;
    const auto f_5 = p / (gamma * q);
    const auto f_6 = 0.5 / gamma;
    const auto f_7 = 0.5 * p / (gamma * q);
    const auto f_8 = gamma / q;
    const auto f_9 = 0.5 / q;
    const auto f_10 = 1.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sph0_0 = buffer.data(sph0 + 0);
    const auto *sph0_3 = buffer.data(sph0 + 3);
    const auto *sph0_5 = buffer.data(sph0 + 5);
    const auto *sph0_6 = buffer.data(sph0 + 6);
    const auto *sph0_9 = buffer.data(sph0 + 9);
    const auto *sph0_24 = buffer.data(sph0 + 24);
    const auto *sph0_27 = buffer.data(sph0 + 27);
    const auto *sph0_36 = buffer.data(sph0 + 36);
    const auto *sph0_38 = buffer.data(sph0 + 38);
    const auto *sph0_39 = buffer.data(sph0 + 39);
    const auto *sph0_41 = buffer.data(sph0 + 41);
    const auto *sph0_42 = buffer.data(sph0 + 42);
    const auto *sph0_47 = buffer.data(sph0 + 47);
    const auto *sph0_51 = buffer.data(sph0 + 51);
    const auto *sph0_57 = buffer.data(sph0 + 57);
    const auto *sph0_59 = buffer.data(sph0 + 59);
    const auto *sph0_60 = buffer.data(sph0 + 60);
    const auto *sph0_62 = buffer.data(sph0 + 62);

    const auto *spg_0 = buffer.data(spg + 0);
    const auto *spg_2 = buffer.data(spg + 2);
    const auto *spg_3 = buffer.data(spg + 3);
    const auto *spg_5 = buffer.data(spg + 5);
    const auto *spg_6 = buffer.data(spg + 6);
    const auto *spg_9 = buffer.data(spg + 9);
    const auto *spg_10 = buffer.data(spg + 10);
    const auto *spg_11 = buffer.data(spg + 11);
    const auto *spg_12 = buffer.data(spg + 12);
    const auto *spg_13 = buffer.data(spg + 13);
    const auto *spg_14 = buffer.data(spg + 14);
    const auto *spg_15 = buffer.data(spg + 15);
    const auto *spg_17 = buffer.data(spg + 17);
    const auto *spg_18 = buffer.data(spg + 18);
    const auto *spg_20 = buffer.data(spg + 20);
    const auto *spg_21 = buffer.data(spg + 21);
    const auto *spg_25 = buffer.data(spg + 25);
    const auto *spg_26 = buffer.data(spg + 26);
    const auto *spg_27 = buffer.data(spg + 27);
    const auto *spg_28 = buffer.data(spg + 28);
    const auto *spg_29 = buffer.data(spg + 29);
    const auto *spg_30 = buffer.data(spg + 30);
    const auto *spg_32 = buffer.data(spg + 32);
    const auto *spg_33 = buffer.data(spg + 33);
    const auto *spg_35 = buffer.data(spg + 35);
    const auto *spg_39 = buffer.data(spg + 39);
    const auto *spg_40 = buffer.data(spg + 40);
    const auto *spg_41 = buffer.data(spg + 41);
    const auto *spg_42 = buffer.data(spg + 42);
    const auto *spg_43 = buffer.data(spg + 43);
    const auto *spg_44 = buffer.data(spg + 44);

    const auto *sph1_0 = buffer.data(sph1 + 0);
    const auto *sph1_3 = buffer.data(sph1 + 3);
    const auto *sph1_5 = buffer.data(sph1 + 5);
    const auto *sph1_6 = buffer.data(sph1 + 6);
    const auto *sph1_9 = buffer.data(sph1 + 9);
    const auto *sph1_24 = buffer.data(sph1 + 24);
    const auto *sph1_27 = buffer.data(sph1 + 27);
    const auto *sph1_36 = buffer.data(sph1 + 36);
    const auto *sph1_38 = buffer.data(sph1 + 38);
    const auto *sph1_39 = buffer.data(sph1 + 39);
    const auto *sph1_41 = buffer.data(sph1 + 41);
    const auto *sph1_42 = buffer.data(sph1 + 42);
    const auto *sph1_47 = buffer.data(sph1 + 47);
    const auto *sph1_51 = buffer.data(sph1 + 51);
    const auto *sph1_57 = buffer.data(sph1 + 57);
    const auto *sph1_59 = buffer.data(sph1 + 59);
    const auto *sph1_60 = buffer.data(sph1 + 60);
    const auto *sph1_62 = buffer.data(sph1 + 62);

    const auto *sdf0_0 = buffer.data(sdf0 + 0);
    const auto *sdf0_3 = buffer.data(sdf0 + 3);
    const auto *sdf0_5 = buffer.data(sdf0 + 5);
    const auto *sdf0_6 = buffer.data(sdf0 + 6);
    const auto *sdf0_8 = buffer.data(sdf0 + 8);
    const auto *sdf0_9 = buffer.data(sdf0 + 9);
    const auto *sdf0_30 = buffer.data(sdf0 + 30);
    const auto *sdf0_33 = buffer.data(sdf0 + 33);
    const auto *sdf0_35 = buffer.data(sdf0 + 35);
    const auto *sdf0_36 = buffer.data(sdf0 + 36);
    const auto *sdf0_38 = buffer.data(sdf0 + 38);
    const auto *sdf0_39 = buffer.data(sdf0 + 39);
    const auto *sdf0_50 = buffer.data(sdf0 + 50);
    const auto *sdf0_53 = buffer.data(sdf0 + 53);
    const auto *sdf0_55 = buffer.data(sdf0 + 55);
    const auto *sdf0_56 = buffer.data(sdf0 + 56);
    const auto *sdf0_58 = buffer.data(sdf0 + 58);
    const auto *sdf0_59 = buffer.data(sdf0 + 59);

    const auto *sdf1_0 = buffer.data(sdf1 + 0);
    const auto *sdf1_3 = buffer.data(sdf1 + 3);
    const auto *sdf1_5 = buffer.data(sdf1 + 5);
    const auto *sdf1_6 = buffer.data(sdf1 + 6);
    const auto *sdf1_8 = buffer.data(sdf1 + 8);
    const auto *sdf1_9 = buffer.data(sdf1 + 9);
    const auto *sdf1_30 = buffer.data(sdf1 + 30);
    const auto *sdf1_33 = buffer.data(sdf1 + 33);
    const auto *sdf1_35 = buffer.data(sdf1 + 35);
    const auto *sdf1_36 = buffer.data(sdf1 + 36);
    const auto *sdf1_38 = buffer.data(sdf1 + 38);
    const auto *sdf1_39 = buffer.data(sdf1 + 39);
    const auto *sdf1_50 = buffer.data(sdf1 + 50);
    const auto *sdf1_53 = buffer.data(sdf1 + 53);
    const auto *sdf1_55 = buffer.data(sdf1 + 55);
    const auto *sdf1_56 = buffer.data(sdf1 + 56);
    const auto *sdf1_58 = buffer.data(sdf1 + 58);
    const auto *sdf1_59 = buffer.data(sdf1 + 59);

    const auto *sdg_0 = buffer.data(sdg + 0);
    const auto *sdg_2 = buffer.data(sdg + 2);
    const auto *sdg_3 = buffer.data(sdg + 3);
    const auto *sdg_5 = buffer.data(sdg + 5);
    const auto *sdg_6 = buffer.data(sdg + 6);
    const auto *sdg_9 = buffer.data(sdg + 9);
    const auto *sdg_10 = buffer.data(sdg + 10);
    const auto *sdg_11 = buffer.data(sdg + 11);
    const auto *sdg_12 = buffer.data(sdg + 12);
    const auto *sdg_13 = buffer.data(sdg + 13);
    const auto *sdg_14 = buffer.data(sdg + 14);
    const auto *sdg_15 = buffer.data(sdg + 15);
    const auto *sdg_17 = buffer.data(sdg + 17);
    const auto *sdg_18 = buffer.data(sdg + 18);
    const auto *sdg_20 = buffer.data(sdg + 20);
    const auto *sdg_25 = buffer.data(sdg + 25);
    const auto *sdg_26 = buffer.data(sdg + 26);
    const auto *sdg_27 = buffer.data(sdg + 27);
    const auto *sdg_28 = buffer.data(sdg + 28);
    const auto *sdg_29 = buffer.data(sdg + 29);
    const auto *sdg_30 = buffer.data(sdg + 30);
    const auto *sdg_32 = buffer.data(sdg + 32);
    const auto *sdg_33 = buffer.data(sdg + 33);
    const auto *sdg_35 = buffer.data(sdg + 35);
    const auto *sdg_40 = buffer.data(sdg + 40);
    const auto *sdg_41 = buffer.data(sdg + 41);
    const auto *sdg_42 = buffer.data(sdg + 42);
    const auto *sdg_43 = buffer.data(sdg + 43);
    const auto *sdg_44 = buffer.data(sdg + 44);
    const auto *sdg_45 = buffer.data(sdg + 45);
    const auto *sdg_47 = buffer.data(sdg + 47);
    const auto *sdg_48 = buffer.data(sdg + 48);
    const auto *sdg_50 = buffer.data(sdg + 50);
    const auto *sdg_51 = buffer.data(sdg + 51);
    const auto *sdg_54 = buffer.data(sdg + 54);
    const auto *sdg_55 = buffer.data(sdg + 55);
    const auto *sdg_56 = buffer.data(sdg + 56);
    const auto *sdg_57 = buffer.data(sdg + 57);
    const auto *sdg_58 = buffer.data(sdg + 58);
    const auto *sdg_59 = buffer.data(sdg + 59);
    const auto *sdg_60 = buffer.data(sdg + 60);
    const auto *sdg_62 = buffer.data(sdg + 62);
    const auto *sdg_63 = buffer.data(sdg + 63);
    const auto *sdg_65 = buffer.data(sdg + 65);
    const auto *sdg_70 = buffer.data(sdg + 70);
    const auto *sdg_71 = buffer.data(sdg + 71);
    const auto *sdg_72 = buffer.data(sdg + 72);
    const auto *sdg_73 = buffer.data(sdg + 73);
    const auto *sdg_74 = buffer.data(sdg + 74);
    const auto *sdg_75 = buffer.data(sdg + 75);
    const auto *sdg_77 = buffer.data(sdg + 77);
    const auto *sdg_78 = buffer.data(sdg + 78);
    const auto *sdg_80 = buffer.data(sdg + 80);
    const auto *sdg_81 = buffer.data(sdg + 81);
    const auto *sdg_84 = buffer.data(sdg + 84);
    const auto *sdg_85 = buffer.data(sdg + 85);
    const auto *sdg_86 = buffer.data(sdg + 86);
    const auto *sdg_87 = buffer.data(sdg + 87);
    const auto *sdg_88 = buffer.data(sdg + 88);
    const auto *sdg_89 = buffer.data(sdg + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pc_x, pc_y, pc_z, spg_0, spg_3, sdf0_0, sdf0_3, \
                         sdf1_0, sdf1_3, sdg_0, sdg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * spg_0[k]
                 + f_1 * sdf0_0[k]
                 - f_2 * sdf1_0[k]
                 + f_3 * pc_x[k] * sdg_0[k];

        t_1[k] = f_3 * pc_y[k] * sdg_0[k];

        t_2[k] = f_3 * pc_z[k] * sdg_0[k];

        t_3[k] = f_0 * spg_3[k]
                 + f_4 * sdf0_3[k]
                 - f_5 * sdf1_3[k]
                 + f_3 * pc_x[k] * sdg_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pc_x, pc_y, spg_5, spg_6, sdf0_5, sdf0_6, sdf1_5, \
                         sdf1_6, sdg_2, sdg_5, sdg_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_3 * pc_y[k] * sdg_2[k];

        t_5[k] = f_0 * spg_5[k]
                 + f_4 * sdf0_5[k]
                 - f_5 * sdf1_5[k]
                 + f_3 * pc_x[k] * sdg_5[k];

        t_6[k] = f_0 * spg_6[k]
                 + f_6 * sdf0_6[k]
                 - f_7 * sdf1_6[k]
                 + f_3 * pc_x[k] * sdg_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, pc_x, pc_y, pc_z, spg_9, spg_10, sdf0_9, sdf1_9, \
                         sdg_3, sdg_5, sdg_9, sdg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_3 * pc_z[k] * sdg_3[k];

        t_8[k] = f_3 * pc_y[k] * sdg_5[k];

        t_9[k] = f_0 * spg_9[k]
                 + f_6 * sdf0_9[k]
                 - f_7 * sdf1_9[k]
                 + f_3 * pc_x[k] * sdg_9[k];

        t_10[k] = f_0 * spg_10[k]
                  + f_3 * pc_x[k] * sdg_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, pc_x, spg_11, spg_12, spg_13, spg_14, sdg_11, \
                         sdg_12, sdg_13, sdg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_0 * spg_11[k]
                  + f_3 * pc_x[k] * sdg_11[k];

        t_12[k] = f_0 * spg_12[k]
                  + f_3 * pc_x[k] * sdg_12[k];

        t_13[k] = f_0 * spg_13[k]
                  + f_3 * pc_x[k] * sdg_13[k];

        t_14[k] = f_0 * spg_14[k]
                  + f_3 * pc_x[k] * sdg_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pc_y, pc_z, sdf0_6, sdf0_8, sdf0_9, sdf1_6, \
                         sdf1_8, sdf1_9, sdg_10, sdg_12, sdg_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * sdf0_6[k]
                  - f_2 * sdf1_6[k]
                  + f_3 * pc_y[k] * sdg_10[k];

        t_16[k] = f_3 * pc_z[k] * sdg_10[k];

        t_17[k] = f_4 * sdf0_8[k]
                  - f_5 * sdf1_8[k]
                  + f_3 * pc_y[k] * sdg_12[k];

        t_18[k] = f_6 * sdf0_9[k]
                  - f_7 * sdf1_9[k]
                  + f_3 * pc_y[k] * sdg_13[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pb_y, pc_y, pc_z, sph0_0, spg_0, \
                         sph1_0, sdf0_9, sdf1_9, sdg_14, sdg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_3 * pc_y[k] * sdg_14[k];

        t_20[k] = f_1 * sdf0_9[k]
                  - f_2 * sdf1_9[k]
                  + f_3 * pc_z[k] * sdg_14[k];

        t_21[k] = pb_y[k] * sph0_0[k]
                  - f_8 * pc_y[k] * sph1_0[k];

        t_22[k] = f_9 * spg_0[k]
                  + f_3 * pc_y[k] * sdg_15[k];

        t_23[k] = f_3 * pc_z[k] * sdg_15[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pb_x, pb_y, pc_x, pc_y, sph0_5, sph0_24, spg_2, \
                         spg_18, sph1_5, sph1_24, sdg_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = pb_x[k] * sph0_24[k]
                  + f_10 * spg_18[k]
                  - f_8 * pc_x[k] * sph1_24[k];

        t_25[k] = f_9 * spg_2[k]
                  + f_3 * pc_y[k] * sdg_17[k];

        t_26[k] = pb_y[k] * sph0_5[k]
                  - f_8 * pc_y[k] * sph1_5[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pb_x, pc_x, pc_y, pc_z, sph0_27, spg_5, spg_21, \
                         sph1_27, sdg_18, sdg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pb_x[k] * sph0_27[k]
                  + f_0 * spg_21[k]
                  - f_8 * pc_x[k] * sph1_27[k];

        t_28[k] = f_3 * pc_z[k] * sdg_18[k];

        t_29[k] = f_9 * spg_5[k]
                  + f_3 * pc_y[k] * sdg_20[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pb_y, pc_x, pc_y, sph0_9, spg_25, spg_26, \
                         spg_27, sph1_9, sdg_25, sdg_26, sdg_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_y[k] * sph0_9[k]
                  - f_8 * pc_y[k] * sph1_9[k];

        t_31[k] = f_9 * spg_25[k]
                  + f_3 * pc_x[k] * sdg_25[k];

        t_32[k] = f_9 * spg_26[k]
                  + f_3 * pc_x[k] * sdg_26[k];

        t_33[k] = f_9 * spg_27[k]
                  + f_3 * pc_x[k] * sdg_27[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, pb_x, pc_x, pc_z, sph0_36, spg_28, spg_29, \
                         sph1_36, sdg_25, sdg_28, sdg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = f_9 * spg_28[k]
                  + f_3 * pc_x[k] * sdg_28[k];

        t_35[k] = f_9 * spg_29[k]
                  + f_3 * pc_x[k] * sdg_29[k];

        t_36[k] = pb_x[k] * sph0_36[k]
                  - f_8 * pc_x[k] * sph1_36[k];

        t_37[k] = f_3 * pc_z[k] * sdg_25[k];
    }

#pragma omp simd aligned(t_38, t_39, t_40, t_41, pb_x, pc_x, pc_y, sph0_38, sph0_39, sph0_41, \
                         spg_14, sph1_38, sph1_39, sph1_41, sdg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_38[k] = pb_x[k] * sph0_38[k]
                  - f_8 * pc_x[k] * sph1_38[k];

        t_39[k] = pb_x[k] * sph0_39[k]
                  - f_8 * pc_x[k] * sph1_39[k];

        t_40[k] = f_9 * spg_14[k]
                  + f_3 * pc_y[k] * sdg_29[k];

        t_41[k] = pb_x[k] * sph0_41[k]
                  - f_8 * pc_x[k] * sph1_41[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pb_z, pc_y, pc_z, sph0_0, sph0_3, \
                         spg_0, sph1_0, sph1_3, sdg_30, sdg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_z[k] * sph0_0[k]
                  - f_8 * pc_z[k] * sph1_0[k];

        t_43[k] = f_3 * pc_y[k] * sdg_30[k];

        t_44[k] = f_9 * spg_0[k]
                  + f_3 * pc_z[k] * sdg_30[k];

        t_45[k] = pb_z[k] * sph0_3[k]
                  - f_8 * pc_z[k] * sph1_3[k];

        t_46[k] = f_3 * pc_y[k] * sdg_32[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, pb_x, pb_z, pc_x, pc_z, sph0_6, sph0_47, spg_3, \
                         spg_35, sph1_6, sph1_47, sdg_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = pb_x[k] * sph0_47[k]
                  + f_10 * spg_35[k]
                  - f_8 * pc_x[k] * sph1_47[k];

        t_48[k] = pb_z[k] * sph0_6[k]
                  - f_8 * pc_z[k] * sph1_6[k];

        t_49[k] = f_9 * spg_3[k]
                  + f_3 * pc_z[k] * sdg_33[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pb_x, pc_x, pc_y, sph0_51, spg_39, spg_40, \
                         spg_41, sph1_51, sdg_35, sdg_40, sdg_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_3 * pc_y[k] * sdg_35[k];

        t_51[k] = pb_x[k] * sph0_51[k]
                  + f_0 * spg_39[k]
                  - f_8 * pc_x[k] * sph1_51[k];

        t_52[k] = f_9 * spg_40[k]
                  + f_3 * pc_x[k] * sdg_40[k];

        t_53[k] = f_9 * spg_41[k]
                  + f_3 * pc_x[k] * sdg_41[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pb_x, pc_x, sph0_57, spg_42, spg_43, spg_44, \
                         sph1_57, sdg_42, sdg_43, sdg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_9 * spg_42[k]
                  + f_3 * pc_x[k] * sdg_42[k];

        t_55[k] = f_9 * spg_43[k]
                  + f_3 * pc_x[k] * sdg_43[k];

        t_56[k] = f_9 * spg_44[k]
                  + f_3 * pc_x[k] * sdg_44[k];

        t_57[k] = pb_x[k] * sph0_57[k]
                  - f_8 * pc_x[k] * sph1_57[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pb_x, pc_x, pc_y, pc_z, sph0_59, sph0_60, \
                         spg_10, sph1_59, sph1_60, sdg_40, sdg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_9 * spg_10[k]
                  + f_3 * pc_z[k] * sdg_40[k];

        t_59[k] = pb_x[k] * sph0_59[k]
                  - f_8 * pc_x[k] * sph1_59[k];

        t_60[k] = pb_x[k] * sph0_60[k]
                  - f_8 * pc_x[k] * sph1_60[k];

        t_61[k] = f_3 * pc_y[k] * sdg_44[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pb_x, pc_x, pc_y, pc_z, sph0_62, spg_15, \
                         sph1_62, sdf0_30, sdf1_30, sdg_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pb_x[k] * sph0_62[k]
                  - f_8 * pc_x[k] * sph1_62[k];

        t_63[k] = f_1 * sdf0_30[k]
                  - f_2 * sdf1_30[k]
                  + f_3 * pc_x[k] * sdg_45[k];

        t_64[k] = f_0 * spg_15[k]
                  + f_3 * pc_y[k] * sdg_45[k];

        t_65[k] = f_3 * pc_z[k] * sdg_45[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pc_x, pc_y, spg_17, sdf0_33, sdf0_35, sdf1_33, \
                         sdf1_35, sdg_47, sdg_48, sdg_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = f_4 * sdf0_33[k]
                  - f_5 * sdf1_33[k]
                  + f_3 * pc_x[k] * sdg_48[k];

        t_67[k] = f_0 * spg_17[k]
                  + f_3 * pc_y[k] * sdg_47[k];

        t_68[k] = f_4 * sdf0_35[k]
                  - f_5 * sdf1_35[k]
                  + f_3 * pc_x[k] * sdg_50[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pc_x, pc_y, pc_z, spg_20, sdf0_36, sdf0_39, \
                         sdf1_36, sdf1_39, sdg_48, sdg_50, sdg_51, \
                         sdg_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = f_6 * sdf0_36[k]
                  - f_7 * sdf1_36[k]
                  + f_3 * pc_x[k] * sdg_51[k];

        t_70[k] = f_3 * pc_z[k] * sdg_48[k];

        t_71[k] = f_0 * spg_20[k]
                  + f_3 * pc_y[k] * sdg_50[k];

        t_72[k] = f_6 * sdf0_39[k]
                  - f_7 * sdf1_39[k]
                  + f_3 * pc_x[k] * sdg_54[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, t_78, pc_x, pc_y, spg_25, sdf0_36, \
                         sdf1_36, sdg_55, sdg_56, sdg_57, sdg_58, \
                         sdg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_3 * pc_x[k] * sdg_55[k];

        t_74[k] = f_3 * pc_x[k] * sdg_56[k];

        t_75[k] = f_3 * pc_x[k] * sdg_57[k];

        t_76[k] = f_3 * pc_x[k] * sdg_58[k];

        t_77[k] = f_3 * pc_x[k] * sdg_59[k];

        t_78[k] = f_0 * spg_25[k]
                  + f_1 * sdf0_36[k]
                  - f_2 * sdf1_36[k]
                  + f_3 * pc_y[k] * sdg_55[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, pc_y, pc_z, spg_27, spg_28, sdf0_38, sdf0_39, \
                         sdf1_38, sdf1_39, sdg_55, sdg_57, sdg_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_3 * pc_z[k] * sdg_55[k];

        t_80[k] = f_0 * spg_27[k]
                  + f_4 * sdf0_38[k]
                  - f_5 * sdf1_38[k]
                  + f_3 * pc_y[k] * sdg_57[k];

        t_81[k] = f_0 * spg_28[k]
                  + f_6 * sdf0_39[k]
                  - f_7 * sdf1_39[k]
                  + f_3 * pc_y[k] * sdg_58[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, t_85, pb_y, pc_y, pc_z, sph0_42, spg_29, spg_30, \
                         sph1_42, sdf0_39, sdf1_39, sdg_59, sdg_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_0 * spg_29[k]
                  + f_3 * pc_y[k] * sdg_59[k];

        t_83[k] = f_1 * sdf0_39[k]
                  - f_2 * sdf1_39[k]
                  + f_3 * pc_z[k] * sdg_59[k];

        t_84[k] = pb_y[k] * sph0_42[k]
                  - f_8 * pc_y[k] * sph1_42[k];

        t_85[k] = f_9 * spg_30[k]
                  + f_3 * pc_y[k] * sdg_60[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, t_89, pb_y, pb_z, pc_y, pc_z, sph0_24, sph0_47, \
                         spg_15, spg_32, sph1_24, sph1_47, sdg_60, \
                         sdg_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_9 * spg_15[k]
                  + f_3 * pc_z[k] * sdg_60[k];

        t_87[k] = pb_z[k] * sph0_24[k]
                  - f_8 * pc_z[k] * sph1_24[k];

        t_88[k] = f_9 * spg_32[k]
                  + f_3 * pc_y[k] * sdg_62[k];

        t_89[k] = pb_y[k] * sph0_47[k]
                  - f_8 * pc_y[k] * sph1_47[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, pb_y, pb_z, pc_y, pc_z, sph0_27, sph0_51, \
                         spg_18, spg_35, sph1_27, sph1_51, sdg_63, \
                         sdg_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = pb_z[k] * sph0_27[k]
                  - f_8 * pc_z[k] * sph1_27[k];

        t_91[k] = f_9 * spg_18[k]
                  + f_3 * pc_z[k] * sdg_63[k];

        t_92[k] = f_9 * spg_35[k]
                  + f_3 * pc_y[k] * sdg_65[k];

        t_93[k] = pb_y[k] * sph0_51[k]
                  - f_8 * pc_y[k] * sph1_51[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, t_97, t_98, t_99, pb_z, pc_x, pc_z, sph0_36, \
                         sph1_36, sdg_70, sdg_71, sdg_72, sdg_73, \
                         sdg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_3 * pc_x[k] * sdg_70[k];

        t_95[k] = f_3 * pc_x[k] * sdg_71[k];

        t_96[k] = f_3 * pc_x[k] * sdg_72[k];

        t_97[k] = f_3 * pc_x[k] * sdg_73[k];

        t_98[k] = f_3 * pc_x[k] * sdg_74[k];

        t_99[k] = pb_z[k] * sph0_36[k]
                  - f_8 * pc_z[k] * sph1_36[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, pb_y, pc_y, pc_z, sph0_59, sph0_60, spg_25, \
                         spg_42, spg_43, sph1_59, sph1_60, sdg_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = f_9 * spg_25[k]
                   + f_3 * pc_z[k] * sdg_70[k];

        t_101[k] = pb_y[k] * sph0_59[k]
                   + f_10 * spg_42[k]
                   - f_8 * pc_y[k] * sph1_59[k];

        t_102[k] = pb_y[k] * sph0_60[k]
                   + f_0 * spg_43[k]
                   - f_8 * pc_y[k] * sph1_60[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, pb_y, pc_x, pc_y, sph0_62, spg_44, \
                         sph1_62, sdf0_50, sdf1_50, sdg_74, sdg_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_9 * spg_44[k]
                   + f_3 * pc_y[k] * sdg_74[k];

        t_104[k] = pb_y[k] * sph0_62[k]
                   - f_8 * pc_y[k] * sph1_62[k];

        t_105[k] = f_1 * sdf0_50[k]
                   - f_2 * sdf1_50[k]
                   + f_3 * pc_x[k] * sdg_75[k];

        t_106[k] = f_3 * pc_y[k] * sdg_75[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pc_x, pc_y, pc_z, spg_30, sdf0_53, \
                         sdf0_55, sdf1_53, sdf1_55, sdg_75, sdg_77, sdg_78, \
                         sdg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_0 * spg_30[k]
                   + f_3 * pc_z[k] * sdg_75[k];

        t_108[k] = f_4 * sdf0_53[k]
                   - f_5 * sdf1_53[k]
                   + f_3 * pc_x[k] * sdg_78[k];

        t_109[k] = f_3 * pc_y[k] * sdg_77[k];

        t_110[k] = f_4 * sdf0_55[k]
                   - f_5 * sdf1_55[k]
                   + f_3 * pc_x[k] * sdg_80[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pc_x, pc_y, pc_z, spg_33, sdf0_56, \
                         sdf0_59, sdf1_56, sdf1_59, sdg_78, sdg_80, sdg_81, \
                         sdg_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_6 * sdf0_56[k]
                   - f_7 * sdf1_56[k]
                   + f_3 * pc_x[k] * sdg_81[k];

        t_112[k] = f_0 * spg_33[k]
                   + f_3 * pc_z[k] * sdg_78[k];

        t_113[k] = f_3 * pc_y[k] * sdg_80[k];

        t_114[k] = f_6 * sdf0_59[k]
                   - f_7 * sdf1_59[k]
                   + f_3 * pc_x[k] * sdg_84[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, t_120, pc_x, pc_y, sdf0_56, \
                         sdf1_56, sdg_85, sdg_86, sdg_87, sdg_88, \
                         sdg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_3 * pc_x[k] * sdg_85[k];

        t_116[k] = f_3 * pc_x[k] * sdg_86[k];

        t_117[k] = f_3 * pc_x[k] * sdg_87[k];

        t_118[k] = f_3 * pc_x[k] * sdg_88[k];

        t_119[k] = f_3 * pc_x[k] * sdg_89[k];

        t_120[k] = f_1 * sdf0_56[k]
                   - f_2 * sdf1_56[k]
                   + f_3 * pc_y[k] * sdg_85[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pc_y, pc_z, spg_40, sdf0_58, sdf0_59, \
                         sdf1_58, sdf1_59, sdg_85, sdg_87, sdg_88, \
                         sdg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_0 * spg_40[k]
                   + f_3 * pc_z[k] * sdg_85[k];

        t_122[k] = f_4 * sdf0_58[k]
                   - f_5 * sdf1_58[k]
                   + f_3 * pc_y[k] * sdg_87[k];

        t_123[k] = f_6 * sdf0_59[k]
                   - f_7 * sdf1_59[k]
                   + f_3 * pc_y[k] * sdg_88[k];

        t_124[k] = f_3 * pc_y[k] * sdg_89[k];
    }

#pragma omp simd aligned(t_125, pc_z, spg_44, sdf0_59, sdf1_59, \
                         sdg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = f_0 * spg_44[k]
                   + f_1 * sdf0_59[k]
                   - f_2 * sdf1_59[k]
                   + f_3 * pc_z[k] * sdg_89[k];
    }
}

}  // namespace simdt3ceri
