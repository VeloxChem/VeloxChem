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


#include "SimdThreeCenterElectronRepulsionVrrRecSDF.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_sdf_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t spf0, const size_t spd,
                                                   const size_t spf1, const size_t sdp0,
                                                   const size_t sdp1, const size_t sdd,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
    const auto f_1 = 1.0 / gamma;
    const auto f_2 = p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = gamma / q;
    const auto f_5 = 0.5 / q;

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

    const auto *spf0_0 = buffer.data(spf0 + 0);
    const auto *spf0_16 = buffer.data(spf0 + 16);
    const auto *spf0_19 = buffer.data(spf0 + 19);
    const auto *spf0_20 = buffer.data(spf0 + 20);
    const auto *spf0_26 = buffer.data(spf0 + 26);
    const auto *spf0_29 = buffer.data(spf0 + 29);

    const auto *spd_0 = buffer.data(spd + 0);
    const auto *spd_3 = buffer.data(spd + 3);
    const auto *spd_4 = buffer.data(spd + 4);
    const auto *spd_5 = buffer.data(spd + 5);
    const auto *spd_6 = buffer.data(spd + 6);
    const auto *spd_9 = buffer.data(spd + 9);
    const auto *spd_10 = buffer.data(spd + 10);
    const auto *spd_11 = buffer.data(spd + 11);
    const auto *spd_12 = buffer.data(spd + 12);
    const auto *spd_15 = buffer.data(spd + 15);
    const auto *spd_16 = buffer.data(spd + 16);
    const auto *spd_17 = buffer.data(spd + 17);

    const auto *spf1_0 = buffer.data(spf1 + 0);
    const auto *spf1_16 = buffer.data(spf1 + 16);
    const auto *spf1_19 = buffer.data(spf1 + 19);
    const auto *spf1_20 = buffer.data(spf1 + 20);
    const auto *spf1_26 = buffer.data(spf1 + 26);
    const auto *spf1_29 = buffer.data(spf1 + 29);

    const auto *sdp0_0 = buffer.data(sdp0 + 0);
    const auto *sdp0_1 = buffer.data(sdp0 + 1);
    const auto *sdp0_2 = buffer.data(sdp0 + 2);
    const auto *sdp0_9 = buffer.data(sdp0 + 9);
    const auto *sdp0_10 = buffer.data(sdp0 + 10);
    const auto *sdp0_11 = buffer.data(sdp0 + 11);
    const auto *sdp0_15 = buffer.data(sdp0 + 15);
    const auto *sdp0_16 = buffer.data(sdp0 + 16);
    const auto *sdp0_17 = buffer.data(sdp0 + 17);

    const auto *sdp1_0 = buffer.data(sdp1 + 0);
    const auto *sdp1_1 = buffer.data(sdp1 + 1);
    const auto *sdp1_2 = buffer.data(sdp1 + 2);
    const auto *sdp1_9 = buffer.data(sdp1 + 9);
    const auto *sdp1_10 = buffer.data(sdp1 + 10);
    const auto *sdp1_11 = buffer.data(sdp1 + 11);
    const auto *sdp1_15 = buffer.data(sdp1 + 15);
    const auto *sdp1_16 = buffer.data(sdp1 + 16);
    const auto *sdp1_17 = buffer.data(sdp1 + 17);

    const auto *sdd_0 = buffer.data(sdd + 0);
    const auto *sdd_3 = buffer.data(sdd + 3);
    const auto *sdd_4 = buffer.data(sdd + 4);
    const auto *sdd_5 = buffer.data(sdd + 5);
    const auto *sdd_6 = buffer.data(sdd + 6);
    const auto *sdd_9 = buffer.data(sdd + 9);
    const auto *sdd_10 = buffer.data(sdd + 10);
    const auto *sdd_11 = buffer.data(sdd + 11);
    const auto *sdd_12 = buffer.data(sdd + 12);
    const auto *sdd_15 = buffer.data(sdd + 15);
    const auto *sdd_16 = buffer.data(sdd + 16);
    const auto *sdd_17 = buffer.data(sdd + 17);
    const auto *sdd_18 = buffer.data(sdd + 18);
    const auto *sdd_21 = buffer.data(sdd + 21);
    const auto *sdd_22 = buffer.data(sdd + 22);
    const auto *sdd_23 = buffer.data(sdd + 23);
    const auto *sdd_24 = buffer.data(sdd + 24);
    const auto *sdd_27 = buffer.data(sdd + 27);
    const auto *sdd_28 = buffer.data(sdd + 28);
    const auto *sdd_29 = buffer.data(sdd + 29);
    const auto *sdd_30 = buffer.data(sdd + 30);
    const auto *sdd_33 = buffer.data(sdd + 33);
    const auto *sdd_34 = buffer.data(sdd + 34);
    const auto *sdd_35 = buffer.data(sdd + 35);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, pc_z, spd_0, spd_3, spd_4, \
                         sdp0_0, sdp1_0, sdd_0, sdd_3, sdd_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * spd_0[k]
                 + f_1 * sdp0_0[k]
                 - f_2 * sdp1_0[k]
                 + f_3 * pc_x[k] * sdd_0[k];

        t_1[k] = f_3 * pc_y[k] * sdd_0[k];

        t_2[k] = f_3 * pc_z[k] * sdd_0[k];

        t_3[k] = f_0 * spd_3[k]
                 + f_3 * pc_x[k] * sdd_3[k];

        t_4[k] = f_0 * spd_4[k]
                 + f_3 * pc_x[k] * sdd_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, pc_x, pc_y, pc_z, spd_5, sdp0_1, sdp0_2, \
                         sdp1_1, sdp1_2, sdd_3, sdd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_0 * spd_5[k]
                 + f_3 * pc_x[k] * sdd_5[k];

        t_6[k] = f_1 * sdp0_1[k]
                 - f_2 * sdp1_1[k]
                 + f_3 * pc_y[k] * sdd_3[k];

        t_7[k] = f_3 * pc_z[k] * sdd_3[k];

        t_8[k] = f_3 * pc_y[k] * sdd_5[k];

        t_9[k] = f_1 * sdp0_2[k]
                 - f_2 * sdp1_2[k]
                 + f_3 * pc_z[k] * sdd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, pb_y, pc_x, pc_y, pc_z, spf0_0, spd_0, spd_9, \
                         spf1_0, sdd_6, sdd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pb_y[k] * spf0_0[k]
                  - f_4 * pc_y[k] * spf1_0[k];

        t_11[k] = f_5 * spd_0[k]
                  + f_3 * pc_y[k] * sdd_6[k];

        t_12[k] = f_3 * pc_z[k] * sdd_6[k];

        t_13[k] = f_5 * spd_9[k]
                  + f_3 * pc_x[k] * sdd_9[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, t_17, pb_x, pc_x, pc_z, spf0_16, spd_10, spd_11, \
                         spf1_16, sdd_9, sdd_10, sdd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_5 * spd_10[k]
                  + f_3 * pc_x[k] * sdd_10[k];

        t_15[k] = f_5 * spd_11[k]
                  + f_3 * pc_x[k] * sdd_11[k];

        t_16[k] = pb_x[k] * spf0_16[k]
                  - f_4 * pc_x[k] * spf1_16[k];

        t_17[k] = f_3 * pc_z[k] * sdd_9[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, pb_x, pb_z, pc_x, pc_y, pc_z, spf0_0, \
                         spf0_19, spd_5, spf1_0, spf1_19, sdd_11, \
                         sdd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * spd_5[k]
                  + f_3 * pc_y[k] * sdd_11[k];

        t_19[k] = pb_x[k] * spf0_19[k]
                  - f_4 * pc_x[k] * spf1_19[k];

        t_20[k] = pb_z[k] * spf0_0[k]
                  - f_4 * pc_z[k] * spf1_0[k];

        t_21[k] = f_3 * pc_y[k] * sdd_12[k];
    }

#pragma omp simd aligned(t_22, t_23, t_24, t_25, pc_x, pc_z, spd_0, spd_15, spd_16, spd_17, \
                         sdd_12, sdd_15, sdd_16, sdd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_22[k] = f_5 * spd_0[k]
                  + f_3 * pc_z[k] * sdd_12[k];

        t_23[k] = f_5 * spd_15[k]
                  + f_3 * pc_x[k] * sdd_15[k];

        t_24[k] = f_5 * spd_16[k]
                  + f_3 * pc_x[k] * sdd_16[k];

        t_25[k] = f_5 * spd_17[k]
                  + f_3 * pc_x[k] * sdd_17[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, t_29, pb_x, pc_x, pc_y, pc_z, spf0_26, spf0_29, \
                         spd_3, spf1_26, spf1_29, sdd_15, sdd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = pb_x[k] * spf0_26[k]
                  - f_4 * pc_x[k] * spf1_26[k];

        t_27[k] = f_5 * spd_3[k]
                  + f_3 * pc_z[k] * sdd_15[k];

        t_28[k] = f_3 * pc_y[k] * sdd_17[k];

        t_29[k] = pb_x[k] * spf0_29[k]
                  - f_4 * pc_x[k] * spf1_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, t_35, pc_x, pc_y, pc_z, spd_6, sdp0_9, \
                         sdp1_9, sdd_18, sdd_21, sdd_22, sdd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_1 * sdp0_9[k]
                  - f_2 * sdp1_9[k]
                  + f_3 * pc_x[k] * sdd_18[k];

        t_31[k] = f_0 * spd_6[k]
                  + f_3 * pc_y[k] * sdd_18[k];

        t_32[k] = f_3 * pc_z[k] * sdd_18[k];

        t_33[k] = f_3 * pc_x[k] * sdd_21[k];

        t_34[k] = f_3 * pc_x[k] * sdd_22[k];

        t_35[k] = f_3 * pc_x[k] * sdd_23[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pc_y, pc_z, spd_9, spd_11, sdp0_10, sdp0_11, \
                         sdp1_10, sdp1_11, sdd_21, sdd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_0 * spd_9[k]
                  + f_1 * sdp0_10[k]
                  - f_2 * sdp1_10[k]
                  + f_3 * pc_y[k] * sdd_21[k];

        t_37[k] = f_3 * pc_z[k] * sdd_21[k];

        t_38[k] = f_0 * spd_11[k]
                  + f_3 * pc_y[k] * sdd_23[k];

        t_39[k] = f_1 * sdp0_11[k]
                  - f_2 * sdp1_11[k]
                  + f_3 * pc_z[k] * sdd_23[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, pb_y, pc_x, pc_y, pc_z, spf0_20, spd_6, \
                         spd_12, spf1_20, sdd_24, sdd_27, sdd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_y[k] * spf0_20[k]
                  - f_4 * pc_y[k] * spf1_20[k];

        t_41[k] = f_5 * spd_12[k]
                  + f_3 * pc_y[k] * sdd_24[k];

        t_42[k] = f_5 * spd_6[k]
                  + f_3 * pc_z[k] * sdd_24[k];

        t_43[k] = f_3 * pc_x[k] * sdd_27[k];

        t_44[k] = f_3 * pc_x[k] * sdd_28[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pb_z, pc_x, pc_y, pc_z, spf0_16, spd_9, \
                         spd_17, spf1_16, sdd_27, sdd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_3 * pc_x[k] * sdd_29[k];

        t_46[k] = pb_z[k] * spf0_16[k]
                  - f_4 * pc_z[k] * spf1_16[k];

        t_47[k] = f_5 * spd_9[k]
                  + f_3 * pc_z[k] * sdd_27[k];

        t_48[k] = f_5 * spd_17[k]
                  + f_3 * pc_y[k] * sdd_29[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, pb_y, pc_x, pc_y, pc_z, spf0_29, \
                         spd_12, spf1_29, sdp0_15, sdp1_15, sdd_30, \
                         sdd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = pb_y[k] * spf0_29[k]
                  - f_4 * pc_y[k] * spf1_29[k];

        t_50[k] = f_1 * sdp0_15[k]
                  - f_2 * sdp1_15[k]
                  + f_3 * pc_x[k] * sdd_30[k];

        t_51[k] = f_3 * pc_y[k] * sdd_30[k];

        t_52[k] = f_0 * spd_12[k]
                  + f_3 * pc_z[k] * sdd_30[k];

        t_53[k] = f_3 * pc_x[k] * sdd_33[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, t_58, pc_x, pc_y, pc_z, spd_15, sdp0_16, \
                         sdp1_16, sdd_33, sdd_34, sdd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_3 * pc_x[k] * sdd_34[k];

        t_55[k] = f_3 * pc_x[k] * sdd_35[k];

        t_56[k] = f_1 * sdp0_16[k]
                  - f_2 * sdp1_16[k]
                  + f_3 * pc_y[k] * sdd_33[k];

        t_57[k] = f_0 * spd_15[k]
                  + f_3 * pc_z[k] * sdd_33[k];

        t_58[k] = f_3 * pc_y[k] * sdd_35[k];
    }

#pragma omp simd aligned(t_59, pc_z, spd_17, sdp0_17, sdp1_17, sdd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_0 * spd_17[k]
                  + f_1 * sdp0_17[k]
                  - f_2 * sdp1_17[k]
                  + f_3 * pc_z[k] * sdd_35[k];
    }
}

}  // namespace simdt3ceri
