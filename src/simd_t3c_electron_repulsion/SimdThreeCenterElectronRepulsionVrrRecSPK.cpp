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


#include "SimdThreeCenterElectronRepulsionVrrRecSPK.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_spk_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t ssk0, const size_t ssi,
                                                   const size_t ssk1, const size_t sph0,
                                                   const size_t sph1, const size_t spi,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / q;
    const auto f_1 = gamma / q;
    const auto f_2 = p / q;
    const auto f_3 = 2.5 / q;
    const auto f_4 = 2.0 / q;
    const auto f_5 = 1.5 / q;
    const auto f_6 = 1.0 / q;
    const auto f_7 = 0.5 / q;
    const auto f_8 = 2.0 / gamma;
    const auto f_9 = 2.0 * p / (gamma * q);
    const auto f_10 = 1.5 / gamma;
    const auto f_11 = 1.5 * p / (gamma * q);
    const auto f_12 = 1.0 / gamma;
    const auto f_13 = p / (gamma * q);
    const auto f_14 = 0.5 / gamma;
    const auto f_15 = 0.5 * p / (gamma * q);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ssk0_0 = buffer.data(ssk0 + 0);
    const auto *ssk0_3 = buffer.data(ssk0 + 3);
    const auto *ssk0_5 = buffer.data(ssk0 + 5);
    const auto *ssk0_6 = buffer.data(ssk0 + 6);
    const auto *ssk0_9 = buffer.data(ssk0 + 9);
    const auto *ssk0_10 = buffer.data(ssk0 + 10);
    const auto *ssk0_12 = buffer.data(ssk0 + 12);
    const auto *ssk0_14 = buffer.data(ssk0 + 14);
    const auto *ssk0_15 = buffer.data(ssk0 + 15);
    const auto *ssk0_17 = buffer.data(ssk0 + 17);
    const auto *ssk0_18 = buffer.data(ssk0 + 18);
    const auto *ssk0_20 = buffer.data(ssk0 + 20);
    const auto *ssk0_28 = buffer.data(ssk0 + 28);
    const auto *ssk0_30 = buffer.data(ssk0 + 30);
    const auto *ssk0_31 = buffer.data(ssk0 + 31);
    const auto *ssk0_32 = buffer.data(ssk0 + 32);
    const auto *ssk0_33 = buffer.data(ssk0 + 33);
    const auto *ssk0_35 = buffer.data(ssk0 + 35);

    const auto *ssi_0 = buffer.data(ssi + 0);
    const auto *ssi_2 = buffer.data(ssi + 2);
    const auto *ssi_3 = buffer.data(ssi + 3);
    const auto *ssi_5 = buffer.data(ssi + 5);
    const auto *ssi_6 = buffer.data(ssi + 6);
    const auto *ssi_9 = buffer.data(ssi + 9);
    const auto *ssi_10 = buffer.data(ssi + 10);
    const auto *ssi_12 = buffer.data(ssi + 12);
    const auto *ssi_14 = buffer.data(ssi + 14);
    const auto *ssi_15 = buffer.data(ssi + 15);
    const auto *ssi_17 = buffer.data(ssi + 17);
    const auto *ssi_18 = buffer.data(ssi + 18);
    const auto *ssi_20 = buffer.data(ssi + 20);
    const auto *ssi_21 = buffer.data(ssi + 21);
    const auto *ssi_22 = buffer.data(ssi + 22);
    const auto *ssi_23 = buffer.data(ssi + 23);
    const auto *ssi_24 = buffer.data(ssi + 24);
    const auto *ssi_25 = buffer.data(ssi + 25);
    const auto *ssi_26 = buffer.data(ssi + 26);
    const auto *ssi_27 = buffer.data(ssi + 27);

    const auto *ssk1_0 = buffer.data(ssk1 + 0);
    const auto *ssk1_3 = buffer.data(ssk1 + 3);
    const auto *ssk1_5 = buffer.data(ssk1 + 5);
    const auto *ssk1_6 = buffer.data(ssk1 + 6);
    const auto *ssk1_9 = buffer.data(ssk1 + 9);
    const auto *ssk1_10 = buffer.data(ssk1 + 10);
    const auto *ssk1_12 = buffer.data(ssk1 + 12);
    const auto *ssk1_14 = buffer.data(ssk1 + 14);
    const auto *ssk1_15 = buffer.data(ssk1 + 15);
    const auto *ssk1_17 = buffer.data(ssk1 + 17);
    const auto *ssk1_18 = buffer.data(ssk1 + 18);
    const auto *ssk1_20 = buffer.data(ssk1 + 20);
    const auto *ssk1_28 = buffer.data(ssk1 + 28);
    const auto *ssk1_30 = buffer.data(ssk1 + 30);
    const auto *ssk1_31 = buffer.data(ssk1 + 31);
    const auto *ssk1_32 = buffer.data(ssk1 + 32);
    const auto *ssk1_33 = buffer.data(ssk1 + 33);
    const auto *ssk1_35 = buffer.data(ssk1 + 35);

    const auto *sph0_24 = buffer.data(sph0 + 24);
    const auto *sph0_27 = buffer.data(sph0 + 27);
    const auto *sph0_31 = buffer.data(sph0 + 31);
    const auto *sph0_33 = buffer.data(sph0 + 33);
    const auto *sph0_36 = buffer.data(sph0 + 36);
    const auto *sph0_38 = buffer.data(sph0 + 38);
    const auto *sph0_39 = buffer.data(sph0 + 39);
    const auto *sph0_47 = buffer.data(sph0 + 47);
    const auto *sph0_51 = buffer.data(sph0 + 51);
    const auto *sph0_54 = buffer.data(sph0 + 54);
    const auto *sph0_56 = buffer.data(sph0 + 56);
    const auto *sph0_59 = buffer.data(sph0 + 59);
    const auto *sph0_60 = buffer.data(sph0 + 60);
    const auto *sph0_61 = buffer.data(sph0 + 61);
    const auto *sph0_62 = buffer.data(sph0 + 62);

    const auto *sph1_24 = buffer.data(sph1 + 24);
    const auto *sph1_27 = buffer.data(sph1 + 27);
    const auto *sph1_31 = buffer.data(sph1 + 31);
    const auto *sph1_33 = buffer.data(sph1 + 33);
    const auto *sph1_36 = buffer.data(sph1 + 36);
    const auto *sph1_38 = buffer.data(sph1 + 38);
    const auto *sph1_39 = buffer.data(sph1 + 39);
    const auto *sph1_47 = buffer.data(sph1 + 47);
    const auto *sph1_51 = buffer.data(sph1 + 51);
    const auto *sph1_54 = buffer.data(sph1 + 54);
    const auto *sph1_56 = buffer.data(sph1 + 56);
    const auto *sph1_59 = buffer.data(sph1 + 59);
    const auto *sph1_60 = buffer.data(sph1 + 60);
    const auto *sph1_61 = buffer.data(sph1 + 61);
    const auto *sph1_62 = buffer.data(sph1 + 62);

    const auto *spi_0 = buffer.data(spi + 0);
    const auto *spi_2 = buffer.data(spi + 2);
    const auto *spi_3 = buffer.data(spi + 3);
    const auto *spi_5 = buffer.data(spi + 5);
    const auto *spi_6 = buffer.data(spi + 6);
    const auto *spi_9 = buffer.data(spi + 9);
    const auto *spi_10 = buffer.data(spi + 10);
    const auto *spi_14 = buffer.data(spi + 14);
    const auto *spi_21 = buffer.data(spi + 21);
    const auto *spi_22 = buffer.data(spi + 22);
    const auto *spi_23 = buffer.data(spi + 23);
    const auto *spi_24 = buffer.data(spi + 24);
    const auto *spi_25 = buffer.data(spi + 25);
    const auto *spi_26 = buffer.data(spi + 26);
    const auto *spi_27 = buffer.data(spi + 27);
    const auto *spi_28 = buffer.data(spi + 28);
    const auto *spi_30 = buffer.data(spi + 30);
    const auto *spi_31 = buffer.data(spi + 31);
    const auto *spi_33 = buffer.data(spi + 33);
    const auto *spi_34 = buffer.data(spi + 34);
    const auto *spi_37 = buffer.data(spi + 37);
    const auto *spi_38 = buffer.data(spi + 38);
    const auto *spi_40 = buffer.data(spi + 40);
    const auto *spi_42 = buffer.data(spi + 42);
    const auto *spi_43 = buffer.data(spi + 43);
    const auto *spi_45 = buffer.data(spi + 45);
    const auto *spi_46 = buffer.data(spi + 46);
    const auto *spi_49 = buffer.data(spi + 49);
    const auto *spi_50 = buffer.data(spi + 50);
    const auto *spi_51 = buffer.data(spi + 51);
    const auto *spi_52 = buffer.data(spi + 52);
    const auto *spi_53 = buffer.data(spi + 53);
    const auto *spi_54 = buffer.data(spi + 54);
    const auto *spi_55 = buffer.data(spi + 55);
    const auto *spi_56 = buffer.data(spi + 56);
    const auto *spi_58 = buffer.data(spi + 58);
    const auto *spi_59 = buffer.data(spi + 59);
    const auto *spi_61 = buffer.data(spi + 61);
    const auto *spi_62 = buffer.data(spi + 62);
    const auto *spi_65 = buffer.data(spi + 65);
    const auto *spi_66 = buffer.data(spi + 66);
    const auto *spi_68 = buffer.data(spi + 68);
    const auto *spi_70 = buffer.data(spi + 70);
    const auto *spi_73 = buffer.data(spi + 73);
    const auto *spi_74 = buffer.data(spi + 74);
    const auto *spi_76 = buffer.data(spi + 76);
    const auto *spi_77 = buffer.data(spi + 77);
    const auto *spi_78 = buffer.data(spi + 78);
    const auto *spi_79 = buffer.data(spi + 79);
    const auto *spi_80 = buffer.data(spi + 80);
    const auto *spi_81 = buffer.data(spi + 81);
    const auto *spi_82 = buffer.data(spi + 82);
    const auto *spi_83 = buffer.data(spi + 83);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pc_x, pc_y, pc_z, ssk0_0, ssk0_3, ssi_0, \
                         ssi_3, ssk1_0, ssk1_3, spi_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = pb_x[k] * ssk0_0[k]
                 + f_0 * ssi_0[k]
                 - f_1 * pc_x[k] * ssk1_0[k];

        t_1[k] = f_2 * pc_y[k] * spi_0[k];

        t_2[k] = f_2 * pc_z[k] * spi_0[k];

        t_3[k] = pb_x[k] * ssk0_3[k]
                 + f_3 * ssi_3[k]
                 - f_1 * pc_x[k] * ssk1_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_x, pc_x, pc_y, pc_z, ssk0_5, ssk0_6, ssi_5, \
                         ssi_6, ssk1_5, ssk1_6, spi_2, spi_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * pc_y[k] * spi_2[k];

        t_5[k] = pb_x[k] * ssk0_5[k]
                 + f_3 * ssi_5[k]
                 - f_1 * pc_x[k] * ssk1_5[k];

        t_6[k] = pb_x[k] * ssk0_6[k]
                 + f_4 * ssi_6[k]
                 - f_1 * pc_x[k] * ssk1_6[k];

        t_7[k] = f_2 * pc_z[k] * spi_3[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_x, pc_x, pc_y, pc_z, ssk0_9, ssk0_10, ssi_9, \
                         ssi_10, ssk1_9, ssk1_10, spi_5, spi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_2 * pc_y[k] * spi_5[k];

        t_9[k] = pb_x[k] * ssk0_9[k]
                 + f_4 * ssi_9[k]
                 - f_1 * pc_x[k] * ssk1_9[k];

        t_10[k] = pb_x[k] * ssk0_10[k]
                  + f_5 * ssi_10[k]
                  - f_1 * pc_x[k] * ssk1_10[k];

        t_11[k] = f_2 * pc_z[k] * spi_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pb_x, pc_x, pc_y, ssk0_12, ssk0_14, ssi_12, ssi_14, \
                         ssk1_12, ssk1_14, spi_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pb_x[k] * ssk0_12[k]
                  + f_5 * ssi_12[k]
                  - f_1 * pc_x[k] * ssk1_12[k];

        t_13[k] = f_2 * pc_y[k] * spi_9[k];

        t_14[k] = pb_x[k] * ssk0_14[k]
                  + f_5 * ssi_14[k]
                  - f_1 * pc_x[k] * ssk1_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, pb_x, pc_x, pc_z, ssk0_15, ssk0_17, ssi_15, ssi_17, \
                         ssk1_15, ssk1_17, spi_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pb_x[k] * ssk0_15[k]
                  + f_6 * ssi_15[k]
                  - f_1 * pc_x[k] * ssk1_15[k];

        t_16[k] = f_2 * pc_z[k] * spi_10[k];

        t_17[k] = pb_x[k] * ssk0_17[k]
                  + f_6 * ssi_17[k]
                  - f_1 * pc_x[k] * ssk1_17[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pb_x, pc_x, pc_y, ssk0_18, ssk0_20, ssi_18, \
                         ssi_20, ssi_21, ssk1_18, ssk1_20, spi_14, \
                         spi_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = pb_x[k] * ssk0_18[k]
                  + f_6 * ssi_18[k]
                  - f_1 * pc_x[k] * ssk1_18[k];

        t_19[k] = f_2 * pc_y[k] * spi_14[k];

        t_20[k] = pb_x[k] * ssk0_20[k]
                  + f_6 * ssi_20[k]
                  - f_1 * pc_x[k] * ssk1_20[k];

        t_21[k] = f_7 * ssi_21[k]
                  + f_2 * pc_x[k] * spi_21[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, t_26, pc_x, ssi_22, ssi_23, ssi_24, ssi_25, \
                         ssi_26, spi_22, spi_23, spi_24, spi_25, \
                         spi_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_7 * ssi_22[k]
                  + f_2 * pc_x[k] * spi_22[k];

        t_23[k] = f_7 * ssi_23[k]
                  + f_2 * pc_x[k] * spi_23[k];

        t_24[k] = f_7 * ssi_24[k]
                  + f_2 * pc_x[k] * spi_24[k];

        t_25[k] = f_7 * ssi_25[k]
                  + f_2 * pc_x[k] * spi_25[k];

        t_26[k] = f_7 * ssi_26[k]
                  + f_2 * pc_x[k] * spi_26[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pb_x, pc_x, pc_z, ssk0_28, ssk0_30, ssi_27, \
                         ssk1_28, ssk1_30, spi_21, spi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_7 * ssi_27[k]
                  + f_2 * pc_x[k] * spi_27[k];

        t_28[k] = pb_x[k] * ssk0_28[k]
                  - f_1 * pc_x[k] * ssk1_28[k];

        t_29[k] = f_2 * pc_z[k] * spi_21[k];

        t_30[k] = pb_x[k] * ssk0_30[k]
                  - f_1 * pc_x[k] * ssk1_30[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, pb_x, pc_x, pc_y, ssk0_31, ssk0_32, ssk0_33, \
                         ssk1_31, ssk1_32, ssk1_33, spi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = pb_x[k] * ssk0_31[k]
                  - f_1 * pc_x[k] * ssk1_31[k];

        t_32[k] = pb_x[k] * ssk0_32[k]
                  - f_1 * pc_x[k] * ssk1_32[k];

        t_33[k] = pb_x[k] * ssk0_33[k]
                  - f_1 * pc_x[k] * ssk1_33[k];

        t_34[k] = f_2 * pc_y[k] * spi_27[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pb_x, pb_y, pc_x, pc_y, pc_z, ssk0_0, \
                         ssk0_35, ssi_0, ssk1_0, ssk1_35, spi_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pb_x[k] * ssk0_35[k]
                  - f_1 * pc_x[k] * ssk1_35[k];

        t_36[k] = pb_y[k] * ssk0_0[k]
                  - f_1 * pc_y[k] * ssk1_0[k];

        t_37[k] = f_7 * ssi_0[k]
                  + f_2 * pc_y[k] * spi_28[k];

        t_38[k] = f_2 * pc_z[k] * spi_28[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_y, pc_x, pc_y, ssk0_5, ssi_2, ssk1_5, sph0_24, \
                         sph1_24, spi_30, spi_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_8 * sph0_24[k]
                  - f_9 * sph1_24[k]
                  + f_2 * pc_x[k] * spi_31[k];

        t_40[k] = f_7 * ssi_2[k]
                  + f_2 * pc_y[k] * spi_30[k];

        t_41[k] = pb_y[k] * ssk0_5[k]
                  - f_1 * pc_y[k] * ssk1_5[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pb_y, pc_x, pc_y, pc_z, ssk0_9, ssi_5, \
                         ssk1_9, sph0_27, sph1_27, spi_31, spi_33, \
                         spi_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_10 * sph0_27[k]
                  - f_11 * sph1_27[k]
                  + f_2 * pc_x[k] * spi_34[k];

        t_43[k] = f_2 * pc_z[k] * spi_31[k];

        t_44[k] = f_7 * ssi_5[k]
                  + f_2 * pc_y[k] * spi_33[k];

        t_45[k] = pb_y[k] * ssk0_9[k]
                  - f_1 * pc_y[k] * ssk1_9[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pc_x, pc_y, pc_z, ssi_9, sph0_31, sph0_33, \
                         sph1_31, sph1_33, spi_34, spi_37, spi_38, \
                         spi_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_12 * sph0_31[k]
                  - f_13 * sph1_31[k]
                  + f_2 * pc_x[k] * spi_38[k];

        t_47[k] = f_2 * pc_z[k] * spi_34[k];

        t_48[k] = f_12 * sph0_33[k]
                  - f_13 * sph1_33[k]
                  + f_2 * pc_x[k] * spi_40[k];

        t_49[k] = f_7 * ssi_9[k]
                  + f_2 * pc_y[k] * spi_37[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pb_y, pc_x, pc_y, pc_z, ssk0_14, ssk1_14, sph0_36, \
                         sph1_36, spi_38, spi_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pb_y[k] * ssk0_14[k]
                  - f_1 * pc_y[k] * ssk1_14[k];

        t_51[k] = f_14 * sph0_36[k]
                  - f_15 * sph1_36[k]
                  + f_2 * pc_x[k] * spi_43[k];

        t_52[k] = f_2 * pc_z[k] * spi_38[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pc_x, pc_y, ssi_14, sph0_38, sph0_39, sph1_38, \
                         sph1_39, spi_42, spi_45, spi_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_14 * sph0_38[k]
                  - f_15 * sph1_38[k]
                  + f_2 * pc_x[k] * spi_45[k];

        t_54[k] = f_14 * sph0_39[k]
                  - f_15 * sph1_39[k]
                  + f_2 * pc_x[k] * spi_46[k];

        t_55[k] = f_7 * ssi_14[k]
                  + f_2 * pc_y[k] * spi_42[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, t_61, pb_y, pc_x, pc_y, ssk0_20, \
                         ssk1_20, spi_49, spi_50, spi_51, spi_52, \
                         spi_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_y[k] * ssk0_20[k]
                  - f_1 * pc_y[k] * ssk1_20[k];

        t_57[k] = f_2 * pc_x[k] * spi_49[k];

        t_58[k] = f_2 * pc_x[k] * spi_50[k];

        t_59[k] = f_2 * pc_x[k] * spi_51[k];

        t_60[k] = f_2 * pc_x[k] * spi_52[k];

        t_61[k] = f_2 * pc_x[k] * spi_53[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pb_y, pc_x, pc_y, pc_z, ssk0_28, ssi_21, \
                         ssk1_28, spi_49, spi_54, spi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = f_2 * pc_x[k] * spi_54[k];

        t_63[k] = f_2 * pc_x[k] * spi_55[k];

        t_64[k] = pb_y[k] * ssk0_28[k]
                  + f_0 * ssi_21[k]
                  - f_1 * pc_y[k] * ssk1_28[k];

        t_65[k] = f_2 * pc_z[k] * spi_49[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, pb_y, pc_y, ssk0_30, ssk0_31, ssk0_32, ssi_23, \
                         ssi_24, ssi_25, ssk1_30, ssk1_31, ssk1_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pb_y[k] * ssk0_30[k]
                  + f_3 * ssi_23[k]
                  - f_1 * pc_y[k] * ssk1_30[k];

        t_67[k] = pb_y[k] * ssk0_31[k]
                  + f_4 * ssi_24[k]
                  - f_1 * pc_y[k] * ssk1_31[k];

        t_68[k] = pb_y[k] * ssk0_32[k]
                  + f_5 * ssi_25[k]
                  - f_1 * pc_y[k] * ssk1_32[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pb_y, pc_y, ssk0_33, ssk0_35, ssi_26, ssi_27, \
                         ssk1_33, ssk1_35, spi_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = pb_y[k] * ssk0_33[k]
                  + f_6 * ssi_26[k]
                  - f_1 * pc_y[k] * ssk1_33[k];

        t_70[k] = f_7 * ssi_27[k]
                  + f_2 * pc_y[k] * spi_55[k];

        t_71[k] = pb_y[k] * ssk0_35[k]
                  - f_1 * pc_y[k] * ssk1_35[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, t_75, t_76, pb_z, pc_y, pc_z, ssk0_0, ssk0_3, \
                         ssi_0, ssk1_0, ssk1_3, spi_56, spi_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = pb_z[k] * ssk0_0[k]
                  - f_1 * pc_z[k] * ssk1_0[k];

        t_73[k] = f_2 * pc_y[k] * spi_56[k];

        t_74[k] = f_7 * ssi_0[k]
                  + f_2 * pc_z[k] * spi_56[k];

        t_75[k] = pb_z[k] * ssk0_3[k]
                  - f_1 * pc_z[k] * ssk1_3[k];

        t_76[k] = f_2 * pc_y[k] * spi_58[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pb_z, pc_x, pc_y, pc_z, ssk0_6, ssi_3, \
                         ssk1_6, sph0_47, sph1_47, spi_59, spi_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_8 * sph0_47[k]
                  - f_9 * sph1_47[k]
                  + f_2 * pc_x[k] * spi_61[k];

        t_78[k] = pb_z[k] * ssk0_6[k]
                  - f_1 * pc_z[k] * ssk1_6[k];

        t_79[k] = f_7 * ssi_3[k]
                  + f_2 * pc_z[k] * spi_59[k];

        t_80[k] = f_2 * pc_y[k] * spi_61[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pb_z, pc_x, pc_z, ssk0_10, ssi_6, ssk1_10, sph0_51, \
                         sph1_51, spi_62, spi_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_10 * sph0_51[k]
                  - f_11 * sph1_51[k]
                  + f_2 * pc_x[k] * spi_65[k];

        t_82[k] = pb_z[k] * ssk0_10[k]
                  - f_1 * pc_z[k] * ssk1_10[k];

        t_83[k] = f_7 * ssi_6[k]
                  + f_2 * pc_z[k] * spi_62[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pc_x, pc_y, sph0_54, sph0_56, sph1_54, sph1_56, \
                         spi_65, spi_68, spi_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_12 * sph0_54[k]
                  - f_13 * sph1_54[k]
                  + f_2 * pc_x[k] * spi_68[k];

        t_85[k] = f_2 * pc_y[k] * spi_65[k];

        t_86[k] = f_12 * sph0_56[k]
                  - f_13 * sph1_56[k]
                  + f_2 * pc_x[k] * spi_70[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pb_z, pc_x, pc_z, ssk0_15, ssi_10, ssk1_15, \
                         sph0_59, sph1_59, spi_66, spi_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = pb_z[k] * ssk0_15[k]
                  - f_1 * pc_z[k] * ssk1_15[k];

        t_88[k] = f_7 * ssi_10[k]
                  + f_2 * pc_z[k] * spi_66[k];

        t_89[k] = f_14 * sph0_59[k]
                  - f_15 * sph1_59[k]
                  + f_2 * pc_x[k] * spi_73[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, t_93, t_94, pc_x, pc_y, sph0_60, sph0_62, sph1_60, \
                         sph1_62, spi_70, spi_74, spi_76, spi_77, \
                         spi_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_14 * sph0_60[k]
                  - f_15 * sph1_60[k]
                  + f_2 * pc_x[k] * spi_74[k];

        t_91[k] = f_2 * pc_y[k] * spi_70[k];

        t_92[k] = f_14 * sph0_62[k]
                  - f_15 * sph1_62[k]
                  + f_2 * pc_x[k] * spi_76[k];

        t_93[k] = f_2 * pc_x[k] * spi_77[k];

        t_94[k] = f_2 * pc_x[k] * spi_78[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, t_99, t_100, pb_z, pc_x, pc_z, ssk0_28, \
                         ssk1_28, spi_79, spi_80, spi_81, spi_82, \
                         spi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_2 * pc_x[k] * spi_79[k];

        t_96[k] = f_2 * pc_x[k] * spi_80[k];

        t_97[k] = f_2 * pc_x[k] * spi_81[k];

        t_98[k] = f_2 * pc_x[k] * spi_82[k];

        t_99[k] = f_2 * pc_x[k] * spi_83[k];

        t_100[k] = pb_z[k] * ssk0_28[k]
                   - f_1 * pc_z[k] * ssk1_28[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, pc_y, pc_z, ssi_21, sph0_59, sph0_60, sph1_59, \
                         sph1_60, spi_77, spi_79, spi_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_7 * ssi_21[k]
                   + f_2 * pc_z[k] * spi_77[k];

        t_102[k] = f_8 * sph0_59[k]
                   - f_9 * sph1_59[k]
                   + f_2 * pc_y[k] * spi_79[k];

        t_103[k] = f_10 * sph0_60[k]
                   - f_11 * sph1_60[k]
                   + f_2 * pc_y[k] * spi_80[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pc_y, sph0_61, sph0_62, sph1_61, sph1_62, \
                         spi_81, spi_82, spi_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = f_12 * sph0_61[k]
                   - f_13 * sph1_61[k]
                   + f_2 * pc_y[k] * spi_81[k];

        t_105[k] = f_14 * sph0_62[k]
                   - f_15 * sph1_62[k]
                   + f_2 * pc_y[k] * spi_82[k];

        t_106[k] = f_2 * pc_y[k] * spi_83[k];
    }

#pragma omp simd aligned(t_107, pb_z, pc_z, ssk0_35, ssi_27, ssk1_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = pb_z[k] * ssk0_35[k]
                   + f_0 * ssi_27[k]
                   - f_1 * pc_z[k] * ssk1_35[k];
    }
}

}  // namespace simdt3ceri
