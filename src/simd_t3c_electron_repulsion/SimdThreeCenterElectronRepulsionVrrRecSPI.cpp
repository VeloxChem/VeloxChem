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


#include "SimdThreeCenterElectronRepulsionVrrRecSPI.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_spi_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t ssi0, const size_t ssh,
                                                   const size_t ssi1, const size_t spg0,
                                                   const size_t spg1, const size_t sph,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / q;
    const auto f_1 = gamma / q;
    const auto f_2 = p / q;
    const auto f_3 = 2.0 / q;
    const auto f_4 = 1.5 / q;
    const auto f_5 = 1.0 / q;
    const auto f_6 = 0.5 / q;
    const auto f_7 = 1.5 / gamma;
    const auto f_8 = 1.5 * p / (gamma * q);
    const auto f_9 = 1.0 / gamma;
    const auto f_10 = p / (gamma * q);
    const auto f_11 = 0.5 / gamma;
    const auto f_12 = 0.5 * p / (gamma * q);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ssi0_0 = buffer.data(ssi0 + 0);
    const auto *ssi0_3 = buffer.data(ssi0 + 3);
    const auto *ssi0_5 = buffer.data(ssi0 + 5);
    const auto *ssi0_6 = buffer.data(ssi0 + 6);
    const auto *ssi0_9 = buffer.data(ssi0 + 9);
    const auto *ssi0_10 = buffer.data(ssi0 + 10);
    const auto *ssi0_12 = buffer.data(ssi0 + 12);
    const auto *ssi0_14 = buffer.data(ssi0 + 14);
    const auto *ssi0_21 = buffer.data(ssi0 + 21);
    const auto *ssi0_23 = buffer.data(ssi0 + 23);
    const auto *ssi0_24 = buffer.data(ssi0 + 24);
    const auto *ssi0_25 = buffer.data(ssi0 + 25);
    const auto *ssi0_27 = buffer.data(ssi0 + 27);

    const auto *ssh_0 = buffer.data(ssh + 0);
    const auto *ssh_2 = buffer.data(ssh + 2);
    const auto *ssh_3 = buffer.data(ssh + 3);
    const auto *ssh_5 = buffer.data(ssh + 5);
    const auto *ssh_6 = buffer.data(ssh + 6);
    const auto *ssh_9 = buffer.data(ssh + 9);
    const auto *ssh_10 = buffer.data(ssh + 10);
    const auto *ssh_12 = buffer.data(ssh + 12);
    const auto *ssh_14 = buffer.data(ssh + 14);
    const auto *ssh_15 = buffer.data(ssh + 15);
    const auto *ssh_16 = buffer.data(ssh + 16);
    const auto *ssh_17 = buffer.data(ssh + 17);
    const auto *ssh_18 = buffer.data(ssh + 18);
    const auto *ssh_19 = buffer.data(ssh + 19);
    const auto *ssh_20 = buffer.data(ssh + 20);

    const auto *ssi1_0 = buffer.data(ssi1 + 0);
    const auto *ssi1_3 = buffer.data(ssi1 + 3);
    const auto *ssi1_5 = buffer.data(ssi1 + 5);
    const auto *ssi1_6 = buffer.data(ssi1 + 6);
    const auto *ssi1_9 = buffer.data(ssi1 + 9);
    const auto *ssi1_10 = buffer.data(ssi1 + 10);
    const auto *ssi1_12 = buffer.data(ssi1 + 12);
    const auto *ssi1_14 = buffer.data(ssi1 + 14);
    const auto *ssi1_21 = buffer.data(ssi1 + 21);
    const auto *ssi1_23 = buffer.data(ssi1 + 23);
    const auto *ssi1_24 = buffer.data(ssi1 + 24);
    const auto *ssi1_25 = buffer.data(ssi1 + 25);
    const auto *ssi1_27 = buffer.data(ssi1 + 27);

    const auto *spg0_18 = buffer.data(spg0 + 18);
    const auto *spg0_21 = buffer.data(spg0 + 21);
    const auto *spg0_25 = buffer.data(spg0 + 25);
    const auto *spg0_27 = buffer.data(spg0 + 27);
    const auto *spg0_35 = buffer.data(spg0 + 35);
    const auto *spg0_39 = buffer.data(spg0 + 39);
    const auto *spg0_42 = buffer.data(spg0 + 42);
    const auto *spg0_43 = buffer.data(spg0 + 43);
    const auto *spg0_44 = buffer.data(spg0 + 44);

    const auto *spg1_18 = buffer.data(spg1 + 18);
    const auto *spg1_21 = buffer.data(spg1 + 21);
    const auto *spg1_25 = buffer.data(spg1 + 25);
    const auto *spg1_27 = buffer.data(spg1 + 27);
    const auto *spg1_35 = buffer.data(spg1 + 35);
    const auto *spg1_39 = buffer.data(spg1 + 39);
    const auto *spg1_42 = buffer.data(spg1 + 42);
    const auto *spg1_43 = buffer.data(spg1 + 43);
    const auto *spg1_44 = buffer.data(spg1 + 44);

    const auto *sph_0 = buffer.data(sph + 0);
    const auto *sph_2 = buffer.data(sph + 2);
    const auto *sph_3 = buffer.data(sph + 3);
    const auto *sph_5 = buffer.data(sph + 5);
    const auto *sph_6 = buffer.data(sph + 6);
    const auto *sph_9 = buffer.data(sph + 9);
    const auto *sph_15 = buffer.data(sph + 15);
    const auto *sph_16 = buffer.data(sph + 16);
    const auto *sph_17 = buffer.data(sph + 17);
    const auto *sph_18 = buffer.data(sph + 18);
    const auto *sph_19 = buffer.data(sph + 19);
    const auto *sph_20 = buffer.data(sph + 20);
    const auto *sph_21 = buffer.data(sph + 21);
    const auto *sph_23 = buffer.data(sph + 23);
    const auto *sph_24 = buffer.data(sph + 24);
    const auto *sph_26 = buffer.data(sph + 26);
    const auto *sph_27 = buffer.data(sph + 27);
    const auto *sph_30 = buffer.data(sph + 30);
    const auto *sph_31 = buffer.data(sph + 31);
    const auto *sph_33 = buffer.data(sph + 33);
    const auto *sph_36 = buffer.data(sph + 36);
    const auto *sph_37 = buffer.data(sph + 37);
    const auto *sph_38 = buffer.data(sph + 38);
    const auto *sph_39 = buffer.data(sph + 39);
    const auto *sph_40 = buffer.data(sph + 40);
    const auto *sph_41 = buffer.data(sph + 41);
    const auto *sph_42 = buffer.data(sph + 42);
    const auto *sph_44 = buffer.data(sph + 44);
    const auto *sph_45 = buffer.data(sph + 45);
    const auto *sph_47 = buffer.data(sph + 47);
    const auto *sph_48 = buffer.data(sph + 48);
    const auto *sph_51 = buffer.data(sph + 51);
    const auto *sph_54 = buffer.data(sph + 54);
    const auto *sph_56 = buffer.data(sph + 56);
    const auto *sph_57 = buffer.data(sph + 57);
    const auto *sph_58 = buffer.data(sph + 58);
    const auto *sph_59 = buffer.data(sph + 59);
    const auto *sph_60 = buffer.data(sph + 60);
    const auto *sph_61 = buffer.data(sph + 61);
    const auto *sph_62 = buffer.data(sph + 62);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pc_x, pc_y, pc_z, ssi0_0, ssi0_3, ssh_0, \
                         ssh_3, ssi1_0, ssi1_3, sph_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = pb_x[k] * ssi0_0[k]
                 + f_0 * ssh_0[k]
                 - f_1 * pc_x[k] * ssi1_0[k];

        t_1[k] = f_2 * pc_y[k] * sph_0[k];

        t_2[k] = f_2 * pc_z[k] * sph_0[k];

        t_3[k] = pb_x[k] * ssi0_3[k]
                 + f_3 * ssh_3[k]
                 - f_1 * pc_x[k] * ssi1_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_x, pc_x, pc_y, pc_z, ssi0_5, ssi0_6, ssh_5, \
                         ssh_6, ssi1_5, ssi1_6, sph_2, sph_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * pc_y[k] * sph_2[k];

        t_5[k] = pb_x[k] * ssi0_5[k]
                 + f_3 * ssh_5[k]
                 - f_1 * pc_x[k] * ssi1_5[k];

        t_6[k] = pb_x[k] * ssi0_6[k]
                 + f_4 * ssh_6[k]
                 - f_1 * pc_x[k] * ssi1_6[k];

        t_7[k] = f_2 * pc_z[k] * sph_3[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_x, pc_x, pc_y, pc_z, ssi0_9, ssi0_10, ssh_9, \
                         ssh_10, ssi1_9, ssi1_10, sph_5, sph_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_2 * pc_y[k] * sph_5[k];

        t_9[k] = pb_x[k] * ssi0_9[k]
                 + f_4 * ssh_9[k]
                 - f_1 * pc_x[k] * ssi1_9[k];

        t_10[k] = pb_x[k] * ssi0_10[k]
                  + f_5 * ssh_10[k]
                  - f_1 * pc_x[k] * ssi1_10[k];

        t_11[k] = f_2 * pc_z[k] * sph_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pb_x, pc_x, pc_y, ssi0_12, ssi0_14, ssh_12, \
                         ssh_14, ssh_15, ssi1_12, ssi1_14, sph_9, \
                         sph_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pb_x[k] * ssi0_12[k]
                  + f_5 * ssh_12[k]
                  - f_1 * pc_x[k] * ssi1_12[k];

        t_13[k] = f_2 * pc_y[k] * sph_9[k];

        t_14[k] = pb_x[k] * ssi0_14[k]
                  + f_5 * ssh_14[k]
                  - f_1 * pc_x[k] * ssi1_14[k];

        t_15[k] = f_6 * ssh_15[k]
                  + f_2 * pc_x[k] * sph_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pc_x, ssh_16, ssh_17, ssh_18, ssh_19, \
                         ssh_20, sph_16, sph_17, sph_18, sph_19, \
                         sph_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_6 * ssh_16[k]
                  + f_2 * pc_x[k] * sph_16[k];

        t_17[k] = f_6 * ssh_17[k]
                  + f_2 * pc_x[k] * sph_17[k];

        t_18[k] = f_6 * ssh_18[k]
                  + f_2 * pc_x[k] * sph_18[k];

        t_19[k] = f_6 * ssh_19[k]
                  + f_2 * pc_x[k] * sph_19[k];

        t_20[k] = f_6 * ssh_20[k]
                  + f_2 * pc_x[k] * sph_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pb_x, pc_x, pc_z, ssi0_21, ssi0_23, ssi0_24, \
                         ssi1_21, ssi1_23, ssi1_24, sph_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = pb_x[k] * ssi0_21[k]
                  - f_1 * pc_x[k] * ssi1_21[k];

        t_22[k] = f_2 * pc_z[k] * sph_15[k];

        t_23[k] = pb_x[k] * ssi0_23[k]
                  - f_1 * pc_x[k] * ssi1_23[k];

        t_24[k] = pb_x[k] * ssi0_24[k]
                  - f_1 * pc_x[k] * ssi1_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pb_x, pb_y, pc_x, pc_y, ssi0_0, ssi0_25, \
                         ssi0_27, ssi1_0, ssi1_25, ssi1_27, sph_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = pb_x[k] * ssi0_25[k]
                  - f_1 * pc_x[k] * ssi1_25[k];

        t_26[k] = f_2 * pc_y[k] * sph_20[k];

        t_27[k] = pb_x[k] * ssi0_27[k]
                  - f_1 * pc_x[k] * ssi1_27[k];

        t_28[k] = pb_y[k] * ssi0_0[k]
                  - f_1 * pc_y[k] * ssi1_0[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pc_x, pc_y, pc_z, ssh_0, ssh_2, spg0_18, \
                         spg1_18, sph_21, sph_23, sph_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_6 * ssh_0[k]
                  + f_2 * pc_y[k] * sph_21[k];

        t_30[k] = f_2 * pc_z[k] * sph_21[k];

        t_31[k] = f_7 * spg0_18[k]
                  - f_8 * spg1_18[k]
                  + f_2 * pc_x[k] * sph_24[k];

        t_32[k] = f_6 * ssh_2[k]
                  + f_2 * pc_y[k] * sph_23[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pb_y, pc_x, pc_y, pc_z, ssi0_5, ssh_5, \
                         ssi1_5, spg0_21, spg1_21, sph_24, sph_26, \
                         sph_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pb_y[k] * ssi0_5[k]
                  - f_1 * pc_y[k] * ssi1_5[k];

        t_34[k] = f_9 * spg0_21[k]
                  - f_10 * spg1_21[k]
                  + f_2 * pc_x[k] * sph_27[k];

        t_35[k] = f_2 * pc_z[k] * sph_24[k];

        t_36[k] = f_6 * ssh_5[k]
                  + f_2 * pc_y[k] * sph_26[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pb_y, pc_x, pc_y, pc_z, ssi0_9, ssi1_9, spg0_25, \
                         spg1_25, sph_27, sph_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = pb_y[k] * ssi0_9[k]
                  - f_1 * pc_y[k] * ssi1_9[k];

        t_38[k] = f_11 * spg0_25[k]
                  - f_12 * spg1_25[k]
                  + f_2 * pc_x[k] * sph_31[k];

        t_39[k] = f_2 * pc_z[k] * sph_27[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pb_y, pc_x, pc_y, ssi0_14, ssh_9, ssi1_14, \
                         spg0_27, spg1_27, sph_30, sph_33, sph_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_11 * spg0_27[k]
                  - f_12 * spg1_27[k]
                  + f_2 * pc_x[k] * sph_33[k];

        t_41[k] = f_6 * ssh_9[k]
                  + f_2 * pc_y[k] * sph_30[k];

        t_42[k] = pb_y[k] * ssi0_14[k]
                  - f_1 * pc_y[k] * ssi1_14[k];

        t_43[k] = f_2 * pc_x[k] * sph_36[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, pc_x, sph_37, sph_38, sph_39, sph_40, \
                         sph_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_2 * pc_x[k] * sph_37[k];

        t_45[k] = f_2 * pc_x[k] * sph_38[k];

        t_46[k] = f_2 * pc_x[k] * sph_39[k];

        t_47[k] = f_2 * pc_x[k] * sph_40[k];

        t_48[k] = f_2 * pc_x[k] * sph_41[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pb_y, pc_y, pc_z, ssi0_21, ssi0_23, ssh_15, ssh_17, \
                         ssi1_21, ssi1_23, sph_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = pb_y[k] * ssi0_21[k]
                  + f_0 * ssh_15[k]
                  - f_1 * pc_y[k] * ssi1_21[k];

        t_50[k] = f_2 * pc_z[k] * sph_36[k];

        t_51[k] = pb_y[k] * ssi0_23[k]
                  + f_3 * ssh_17[k]
                  - f_1 * pc_y[k] * ssi1_23[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, t_55, pb_y, pc_y, ssi0_24, ssi0_25, ssi0_27, \
                         ssh_18, ssh_19, ssh_20, ssi1_24, ssi1_25, ssi1_27, \
                         sph_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = pb_y[k] * ssi0_24[k]
                  + f_4 * ssh_18[k]
                  - f_1 * pc_y[k] * ssi1_24[k];

        t_53[k] = pb_y[k] * ssi0_25[k]
                  + f_5 * ssh_19[k]
                  - f_1 * pc_y[k] * ssi1_25[k];

        t_54[k] = f_6 * ssh_20[k]
                  + f_2 * pc_y[k] * sph_41[k];

        t_55[k] = pb_y[k] * ssi0_27[k]
                  - f_1 * pc_y[k] * ssi1_27[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, pb_z, pc_y, pc_z, ssi0_0, ssi0_3, \
                         ssh_0, ssi1_0, ssi1_3, sph_42, sph_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pb_z[k] * ssi0_0[k]
                  - f_1 * pc_z[k] * ssi1_0[k];

        t_57[k] = f_2 * pc_y[k] * sph_42[k];

        t_58[k] = f_6 * ssh_0[k]
                  + f_2 * pc_z[k] * sph_42[k];

        t_59[k] = pb_z[k] * ssi0_3[k]
                  - f_1 * pc_z[k] * ssi1_3[k];

        t_60[k] = f_2 * pc_y[k] * sph_44[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pb_z, pc_x, pc_y, pc_z, ssi0_6, ssh_3, \
                         ssi1_6, spg0_35, spg1_35, sph_45, sph_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_7 * spg0_35[k]
                  - f_8 * spg1_35[k]
                  + f_2 * pc_x[k] * sph_47[k];

        t_62[k] = pb_z[k] * ssi0_6[k]
                  - f_1 * pc_z[k] * ssi1_6[k];

        t_63[k] = f_6 * ssh_3[k]
                  + f_2 * pc_z[k] * sph_45[k];

        t_64[k] = f_2 * pc_y[k] * sph_47[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, pb_z, pc_x, pc_z, ssi0_10, ssh_6, ssi1_10, spg0_39, \
                         spg1_39, sph_48, sph_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_9 * spg0_39[k]
                  - f_10 * spg1_39[k]
                  + f_2 * pc_x[k] * sph_51[k];

        t_66[k] = pb_z[k] * ssi0_10[k]
                  - f_1 * pc_z[k] * ssi1_10[k];

        t_67[k] = f_6 * ssh_6[k]
                  + f_2 * pc_z[k] * sph_48[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pc_x, pc_y, spg0_42, spg0_44, spg1_42, \
                         spg1_44, sph_51, sph_54, sph_56, sph_57, \
                         sph_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_11 * spg0_42[k]
                  - f_12 * spg1_42[k]
                  + f_2 * pc_x[k] * sph_54[k];

        t_69[k] = f_2 * pc_y[k] * sph_51[k];

        t_70[k] = f_11 * spg0_44[k]
                  - f_12 * spg1_44[k]
                  + f_2 * pc_x[k] * sph_56[k];

        t_71[k] = f_2 * pc_x[k] * sph_57[k];

        t_72[k] = f_2 * pc_x[k] * sph_58[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pb_z, pc_x, pc_z, ssi0_21, ssi1_21, \
                         sph_59, sph_60, sph_61, sph_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_2 * pc_x[k] * sph_59[k];

        t_74[k] = f_2 * pc_x[k] * sph_60[k];

        t_75[k] = f_2 * pc_x[k] * sph_61[k];

        t_76[k] = f_2 * pc_x[k] * sph_62[k];

        t_77[k] = pb_z[k] * ssi0_21[k]
                  - f_1 * pc_z[k] * ssi1_21[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pc_y, pc_z, ssh_15, spg0_42, spg0_43, spg1_42, \
                         spg1_43, sph_57, sph_59, sph_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_6 * ssh_15[k]
                  + f_2 * pc_z[k] * sph_57[k];

        t_79[k] = f_7 * spg0_42[k]
                  - f_8 * spg1_42[k]
                  + f_2 * pc_y[k] * sph_59[k];

        t_80[k] = f_9 * spg0_43[k]
                  - f_10 * spg1_43[k]
                  + f_2 * pc_y[k] * sph_60[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pb_z, pc_y, pc_z, ssi0_27, ssh_20, ssi1_27, \
                         spg0_44, spg1_44, sph_61, sph_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_11 * spg0_44[k]
                  - f_12 * spg1_44[k]
                  + f_2 * pc_y[k] * sph_61[k];

        t_82[k] = f_2 * pc_y[k] * sph_62[k];

        t_83[k] = pb_z[k] * ssi0_27[k]
                  + f_0 * ssh_20[k]
                  - f_1 * pc_z[k] * ssi1_27[k];
    }
}

}  // namespace simdt3ceri
