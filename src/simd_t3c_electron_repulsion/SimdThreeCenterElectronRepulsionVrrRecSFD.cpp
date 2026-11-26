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


#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_sfd_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t sdd0, const size_t sdp,
                                                   const size_t sdd1, const size_t sfs0,
                                                   const size_t sfs1, const size_t sfp,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
    const auto f_1 = 0.5 / gamma;
    const auto f_2 = 0.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 1.0 / q;
    const auto f_6 = 0.5 / q;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sdd0_0 = buffer.data(sdd0 + 0);
    const auto *sdd0_3 = buffer.data(sdd0 + 3);
    const auto *sdd0_5 = buffer.data(sdd0 + 5);
    const auto *sdd0_12 = buffer.data(sdd0 + 12);
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
    const auto *sdp_4 = buffer.data(sdp + 4);
    const auto *sdp_5 = buffer.data(sdp + 5);
    const auto *sdp_7 = buffer.data(sdp + 7);
    const auto *sdp_8 = buffer.data(sdp + 8);
    const auto *sdp_9 = buffer.data(sdp + 9);
    const auto *sdp_10 = buffer.data(sdp + 10);
    const auto *sdp_11 = buffer.data(sdp + 11);
    const auto *sdp_13 = buffer.data(sdp + 13);
    const auto *sdp_14 = buffer.data(sdp + 14);
    const auto *sdp_15 = buffer.data(sdp + 15);
    const auto *sdp_16 = buffer.data(sdp + 16);
    const auto *sdp_17 = buffer.data(sdp + 17);

    const auto *sdd1_0 = buffer.data(sdd1 + 0);
    const auto *sdd1_3 = buffer.data(sdd1 + 3);
    const auto *sdd1_5 = buffer.data(sdd1 + 5);
    const auto *sdd1_12 = buffer.data(sdd1 + 12);
    const auto *sdd1_18 = buffer.data(sdd1 + 18);
    const auto *sdd1_21 = buffer.data(sdd1 + 21);
    const auto *sdd1_23 = buffer.data(sdd1 + 23);
    const auto *sdd1_27 = buffer.data(sdd1 + 27);
    const auto *sdd1_29 = buffer.data(sdd1 + 29);
    const auto *sdd1_30 = buffer.data(sdd1 + 30);
    const auto *sdd1_33 = buffer.data(sdd1 + 33);
    const auto *sdd1_35 = buffer.data(sdd1 + 35);

    const auto *sfs0_0 = buffer.data(sfs0 + 0);
    const auto *sfs0_1 = buffer.data(sfs0 + 1);
    const auto *sfs0_2 = buffer.data(sfs0 + 2);
    const auto *sfs0_6 = buffer.data(sfs0 + 6);
    const auto *sfs0_7 = buffer.data(sfs0 + 7);
    const auto *sfs0_9 = buffer.data(sfs0 + 9);

    const auto *sfs1_0 = buffer.data(sfs1 + 0);
    const auto *sfs1_1 = buffer.data(sfs1 + 1);
    const auto *sfs1_2 = buffer.data(sfs1 + 2);
    const auto *sfs1_6 = buffer.data(sfs1 + 6);
    const auto *sfs1_7 = buffer.data(sfs1 + 7);
    const auto *sfs1_9 = buffer.data(sfs1 + 9);

    const auto *sfp_0 = buffer.data(sfp + 0);
    const auto *sfp_1 = buffer.data(sfp + 1);
    const auto *sfp_2 = buffer.data(sfp + 2);
    const auto *sfp_4 = buffer.data(sfp + 4);
    const auto *sfp_5 = buffer.data(sfp + 5);
    const auto *sfp_7 = buffer.data(sfp + 7);
    const auto *sfp_8 = buffer.data(sfp + 8);
    const auto *sfp_10 = buffer.data(sfp + 10);
    const auto *sfp_11 = buffer.data(sfp + 11);
    const auto *sfp_13 = buffer.data(sfp + 13);
    const auto *sfp_14 = buffer.data(sfp + 14);
    const auto *sfp_16 = buffer.data(sfp + 16);
    const auto *sfp_17 = buffer.data(sfp + 17);
    const auto *sfp_18 = buffer.data(sfp + 18);
    const auto *sfp_19 = buffer.data(sfp + 19);
    const auto *sfp_20 = buffer.data(sfp + 20);
    const auto *sfp_22 = buffer.data(sfp + 22);
    const auto *sfp_23 = buffer.data(sfp + 23);
    const auto *sfp_25 = buffer.data(sfp + 25);
    const auto *sfp_26 = buffer.data(sfp + 26);
    const auto *sfp_27 = buffer.data(sfp + 27);
    const auto *sfp_28 = buffer.data(sfp + 28);
    const auto *sfp_29 = buffer.data(sfp + 29);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, sdp_0, sdp_1, sdp_2, sfs0_0, \
                         sfs1_0, sfp_0, sfp_1, sfp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sdp_0[k]
                 + f_1 * sfs0_0[k]
                 - f_2 * sfs1_0[k]
                 + f_3 * pc_x[k] * sfp_0[k];

        t_1[k] = f_0 * sdp_1[k]
                 + f_3 * pc_x[k] * sfp_1[k];

        t_2[k] = f_0 * sdp_2[k]
                 + f_3 * pc_x[k] * sfp_2[k];

        t_3[k] = f_1 * sfs0_0[k]
                 - f_2 * sfs1_0[k]
                 + f_3 * pc_y[k] * sfp_1[k];

        t_4[k] = f_3 * pc_y[k] * sfp_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, pb_y, pc_x, pc_y, pc_z, sdd0_0, sdp_4, sdd1_0, sfs0_0, \
                         sfs1_0, sfp_2, sfp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_1 * sfs0_0[k]
                 - f_2 * sfs1_0[k]
                 + f_3 * pc_z[k] * sfp_2[k];

        t_6[k] = pb_y[k] * sdd0_0[k]
                 - f_4 * pc_y[k] * sdd1_0[k];

        t_7[k] = f_5 * sdp_4[k]
                 + f_3 * pc_x[k] * sfp_4[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_y, pc_x, pc_y, sdd0_5, sdp_1, sdp_2, sdp_5, \
                         sdd1_5, sfs0_1, sfs1_1, sfp_4, sfp_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_5 * sdp_5[k]
                 + f_3 * pc_x[k] * sfp_5[k];

        t_9[k] = f_6 * sdp_1[k]
                 + f_1 * sfs0_1[k]
                 - f_2 * sfs1_1[k]
                 + f_3 * pc_y[k] * sfp_4[k];

        t_10[k] = f_6 * sdp_2[k]
                  + f_3 * pc_y[k] * sfp_5[k];

        t_11[k] = pb_y[k] * sdd0_5[k]
                  - f_4 * pc_y[k] * sdd1_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pb_z, pc_x, pc_z, sdd0_0, sdd0_3, sdp_7, \
                         sdp_8, sdd1_0, sdd1_3, sfp_7, sfp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pb_z[k] * sdd0_0[k]
                  - f_4 * pc_z[k] * sdd1_0[k];

        t_13[k] = f_5 * sdp_7[k]
                  + f_3 * pc_x[k] * sfp_7[k];

        t_14[k] = f_5 * sdp_8[k]
                  + f_3 * pc_x[k] * sfp_8[k];

        t_15[k] = pb_z[k] * sdd0_3[k]
                  - f_4 * pc_z[k] * sdd1_3[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pb_x, pc_x, pc_y, pc_z, sdd0_18, sdp_2, sdp_9, \
                         sdd1_18, sfs0_2, sfs1_2, sfp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_3 * pc_y[k] * sfp_8[k];

        t_17[k] = f_6 * sdp_2[k]
                  + f_1 * sfs0_2[k]
                  - f_2 * sfs1_2[k]
                  + f_3 * pc_z[k] * sfp_8[k];

        t_18[k] = pb_x[k] * sdd0_18[k]
                  + f_5 * sdp_9[k]
                  - f_4 * pc_x[k] * sdd1_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pb_x, pc_x, pc_y, sdd0_21, sdp_5, sdp_10, \
                         sdp_11, sdd1_21, sfp_10, sfp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_6 * sdp_10[k]
                  + f_3 * pc_x[k] * sfp_10[k];

        t_20[k] = f_6 * sdp_11[k]
                  + f_3 * pc_x[k] * sfp_11[k];

        t_21[k] = pb_x[k] * sdd0_21[k]
                  - f_4 * pc_x[k] * sdd1_21[k];

        t_22[k] = f_5 * sdp_5[k]
                  + f_3 * pc_y[k] * sfp_11[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pb_x, pb_y, pc_x, pc_y, sdd0_12, sdd0_23, \
                         sdp_13, sdp_14, sdd1_12, sdd1_23, sfp_13, \
                         sfp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = pb_x[k] * sdd0_23[k]
                  - f_4 * pc_x[k] * sdd1_23[k];

        t_24[k] = pb_y[k] * sdd0_12[k]
                  - f_4 * pc_y[k] * sdd1_12[k];

        t_25[k] = f_6 * sdp_13[k]
                  + f_3 * pc_x[k] * sfp_13[k];

        t_26[k] = f_6 * sdp_14[k]
                  + f_3 * pc_x[k] * sfp_14[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pb_x, pc_x, pc_y, sdd0_27, sdd0_29, sdd0_30, \
                         sdp_8, sdp_15, sdd1_27, sdd1_29, sdd1_30, \
                         sfp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pb_x[k] * sdd0_27[k]
                  - f_4 * pc_x[k] * sdd1_27[k];

        t_28[k] = f_6 * sdp_8[k]
                  + f_3 * pc_y[k] * sfp_14[k];

        t_29[k] = pb_x[k] * sdd0_29[k]
                  - f_4 * pc_x[k] * sdd1_29[k];

        t_30[k] = pb_x[k] * sdd0_30[k]
                  + f_5 * sdp_15[k]
                  - f_4 * pc_x[k] * sdd1_30[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pb_x, pc_x, pc_y, sdd0_33, sdd0_35, \
                         sdp_16, sdp_17, sdd1_33, sdd1_35, sfp_16, \
                         sfp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_6 * sdp_16[k]
                  + f_3 * pc_x[k] * sfp_16[k];

        t_32[k] = f_6 * sdp_17[k]
                  + f_3 * pc_x[k] * sfp_17[k];

        t_33[k] = pb_x[k] * sdd0_33[k]
                  - f_4 * pc_x[k] * sdd1_33[k];

        t_34[k] = f_3 * pc_y[k] * sfp_17[k];

        t_35[k] = pb_x[k] * sdd0_35[k]
                  - f_4 * pc_x[k] * sdd1_35[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, t_40, t_41, pc_x, pc_y, pc_z, sdp_10, sdp_11, \
                         sfs0_6, sfs1_6, sfp_18, sfp_19, sfp_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_1 * sfs0_6[k]
                  - f_2 * sfs1_6[k]
                  + f_3 * pc_x[k] * sfp_18[k];

        t_37[k] = f_3 * pc_x[k] * sfp_19[k];

        t_38[k] = f_3 * pc_x[k] * sfp_20[k];

        t_39[k] = f_0 * sdp_10[k]
                  + f_1 * sfs0_6[k]
                  - f_2 * sfs1_6[k]
                  + f_3 * pc_y[k] * sfp_19[k];

        t_40[k] = f_0 * sdp_11[k]
                  + f_3 * pc_y[k] * sfp_20[k];

        t_41[k] = f_1 * sfs0_6[k]
                  - f_2 * sfs1_6[k]
                  + f_3 * pc_z[k] * sfp_20[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pb_z, pc_x, pc_y, pc_z, sdd0_18, \
                         sdd0_21, sdp_14, sdd1_18, sdd1_21, sfp_22, \
                         sfp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_z[k] * sdd0_18[k]
                  - f_4 * pc_z[k] * sdd1_18[k];

        t_43[k] = f_3 * pc_x[k] * sfp_22[k];

        t_44[k] = f_3 * pc_x[k] * sfp_23[k];

        t_45[k] = pb_z[k] * sdd0_21[k]
                  - f_4 * pc_z[k] * sdd1_21[k];

        t_46[k] = f_5 * sdp_14[k]
                  + f_3 * pc_y[k] * sfp_23[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pb_y, pc_x, pc_y, pc_z, sdd0_30, sdp_11, \
                         sdd1_30, sfs0_7, sfs1_7, sfp_23, sfp_25, \
                         sfp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_6 * sdp_11[k]
                  + f_1 * sfs0_7[k]
                  - f_2 * sfs1_7[k]
                  + f_3 * pc_z[k] * sfp_23[k];

        t_48[k] = pb_y[k] * sdd0_30[k]
                  - f_4 * pc_y[k] * sdd1_30[k];

        t_49[k] = f_3 * pc_x[k] * sfp_25[k];

        t_50[k] = f_3 * pc_x[k] * sfp_26[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pb_y, pc_y, sdd0_33, sdd0_35, sdp_16, sdp_17, \
                         sdd1_33, sdd1_35, sfp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pb_y[k] * sdd0_33[k]
                  + f_5 * sdp_16[k]
                  - f_4 * pc_y[k] * sdd1_33[k];

        t_52[k] = f_6 * sdp_17[k]
                  + f_3 * pc_y[k] * sfp_26[k];

        t_53[k] = pb_y[k] * sdd0_35[k]
                  - f_4 * pc_y[k] * sdd1_35[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, t_59, pc_x, pc_y, pc_z, sdp_17, sfs0_9, \
                         sfs1_9, sfp_27, sfp_28, sfp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_1 * sfs0_9[k]
                  - f_2 * sfs1_9[k]
                  + f_3 * pc_x[k] * sfp_27[k];

        t_55[k] = f_3 * pc_x[k] * sfp_28[k];

        t_56[k] = f_3 * pc_x[k] * sfp_29[k];

        t_57[k] = f_1 * sfs0_9[k]
                  - f_2 * sfs1_9[k]
                  + f_3 * pc_y[k] * sfp_28[k];

        t_58[k] = f_3 * pc_y[k] * sfp_29[k];

        t_59[k] = f_0 * sdp_17[k]
                  + f_1 * sfs0_9[k]
                  - f_2 * sfs1_9[k]
                  + f_3 * pc_z[k] * sfp_29[k];
    }
}

}  // namespace simdt3ceri
