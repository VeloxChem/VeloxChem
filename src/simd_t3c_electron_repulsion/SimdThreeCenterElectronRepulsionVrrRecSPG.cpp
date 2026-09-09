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


#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_spg_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t ssg0, const size_t ssf,
                                                   const size_t ssg1, const size_t spd0,
                                                   const size_t spd1, const size_t spf,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / q;
    const auto f_1 = gamma / q;
    const auto f_2 = p / q;
    const auto f_3 = 1.0 / q;
    const auto f_4 = 0.5 / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ssg0_0 = buffer.data(ssg0 + 0);
    const auto *ssg0_3 = buffer.data(ssg0 + 3);
    const auto *ssg0_5 = buffer.data(ssg0 + 5);
    const auto *ssg0_10 = buffer.data(ssg0 + 10);
    const auto *ssg0_12 = buffer.data(ssg0 + 12);
    const auto *ssg0_14 = buffer.data(ssg0 + 14);

    const auto *ssf_0 = buffer.data(ssf + 0);
    const auto *ssf_2 = buffer.data(ssf + 2);
    const auto *ssf_3 = buffer.data(ssf + 3);
    const auto *ssf_5 = buffer.data(ssf + 5);
    const auto *ssf_6 = buffer.data(ssf + 6);
    const auto *ssf_7 = buffer.data(ssf + 7);
    const auto *ssf_8 = buffer.data(ssf + 8);
    const auto *ssf_9 = buffer.data(ssf + 9);

    const auto *ssg1_0 = buffer.data(ssg1 + 0);
    const auto *ssg1_3 = buffer.data(ssg1 + 3);
    const auto *ssg1_5 = buffer.data(ssg1 + 5);
    const auto *ssg1_10 = buffer.data(ssg1 + 10);
    const auto *ssg1_12 = buffer.data(ssg1 + 12);
    const auto *ssg1_14 = buffer.data(ssg1 + 14);

    const auto *spd0_9 = buffer.data(spd0 + 9);
    const auto *spd0_17 = buffer.data(spd0 + 17);

    const auto *spd1_9 = buffer.data(spd1 + 9);
    const auto *spd1_17 = buffer.data(spd1 + 17);

    const auto *spf_0 = buffer.data(spf + 0);
    const auto *spf_2 = buffer.data(spf + 2);
    const auto *spf_6 = buffer.data(spf + 6);
    const auto *spf_7 = buffer.data(spf + 7);
    const auto *spf_8 = buffer.data(spf + 8);
    const auto *spf_9 = buffer.data(spf + 9);
    const auto *spf_10 = buffer.data(spf + 10);
    const auto *spf_12 = buffer.data(spf + 12);
    const auto *spf_13 = buffer.data(spf + 13);
    const auto *spf_16 = buffer.data(spf + 16);
    const auto *spf_17 = buffer.data(spf + 17);
    const auto *spf_18 = buffer.data(spf + 18);
    const auto *spf_19 = buffer.data(spf + 19);
    const auto *spf_20 = buffer.data(spf + 20);
    const auto *spf_22 = buffer.data(spf + 22);
    const auto *spf_25 = buffer.data(spf + 25);
    const auto *spf_26 = buffer.data(spf + 26);
    const auto *spf_27 = buffer.data(spf + 27);
    const auto *spf_28 = buffer.data(spf + 28);
    const auto *spf_29 = buffer.data(spf + 29);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pc_x, pc_y, pc_z, ssg0_0, ssg0_3, ssf_0, \
                         ssf_3, ssg1_0, ssg1_3, spf_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = pb_x[k] * ssg0_0[k]
                 + f_0 * ssf_0[k]
                 - f_1 * pc_x[k] * ssg1_0[k];

        t_1[k] = f_2 * pc_y[k] * spf_0[k];

        t_2[k] = f_2 * pc_z[k] * spf_0[k];

        t_3[k] = pb_x[k] * ssg0_3[k]
                 + f_3 * ssf_3[k]
                 - f_1 * pc_x[k] * ssg1_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_x, pc_x, pc_y, ssg0_5, ssf_5, ssf_6, ssf_7, \
                         ssg1_5, spf_2, spf_6, spf_7 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * pc_y[k] * spf_2[k];

        t_5[k] = pb_x[k] * ssg0_5[k]
                 + f_3 * ssf_5[k]
                 - f_1 * pc_x[k] * ssg1_5[k];

        t_6[k] = f_4 * ssf_6[k]
                 + f_2 * pc_x[k] * spf_6[k];

        t_7[k] = f_4 * ssf_7[k]
                 + f_2 * pc_x[k] * spf_7[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_x, pc_x, pc_z, ssg0_10, ssf_8, ssf_9, \
                         ssg1_10, spf_6, spf_8, spf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_4 * ssf_8[k]
                 + f_2 * pc_x[k] * spf_8[k];

        t_9[k] = f_4 * ssf_9[k]
                 + f_2 * pc_x[k] * spf_9[k];

        t_10[k] = pb_x[k] * ssg0_10[k]
                  - f_1 * pc_x[k] * ssg1_10[k];

        t_11[k] = f_2 * pc_z[k] * spf_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pb_x, pb_y, pc_x, pc_y, ssg0_0, ssg0_12, \
                         ssg0_14, ssg1_0, ssg1_12, ssg1_14, spf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pb_x[k] * ssg0_12[k]
                  - f_1 * pc_x[k] * ssg1_12[k];

        t_13[k] = f_2 * pc_y[k] * spf_9[k];

        t_14[k] = pb_x[k] * ssg0_14[k]
                  - f_1 * pc_x[k] * ssg1_14[k];

        t_15[k] = pb_y[k] * ssg0_0[k]
                  - f_1 * pc_y[k] * ssg1_0[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pc_x, pc_y, pc_z, ssf_0, ssf_2, spd0_9, \
                         spd1_9, spf_10, spf_12, spf_13 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_4 * ssf_0[k]
                  + f_2 * pc_y[k] * spf_10[k];

        t_17[k] = f_2 * pc_z[k] * spf_10[k];

        t_18[k] = f_5 * spd0_9[k]
                  - f_6 * spd1_9[k]
                  + f_2 * pc_x[k] * spf_13[k];

        t_19[k] = f_4 * ssf_2[k]
                  + f_2 * pc_y[k] * spf_12[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, pb_y, pc_x, pc_y, ssg0_5, ssg1_5, \
                         spf_16, spf_17, spf_18, spf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pb_y[k] * ssg0_5[k]
                  - f_1 * pc_y[k] * ssg1_5[k];

        t_21[k] = f_2 * pc_x[k] * spf_16[k];

        t_22[k] = f_2 * pc_x[k] * spf_17[k];

        t_23[k] = f_2 * pc_x[k] * spf_18[k];

        t_24[k] = f_2 * pc_x[k] * spf_19[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, pb_y, pc_y, pc_z, ssg0_10, ssg0_12, ssf_6, \
                         ssf_8, ssf_9, ssg1_10, ssg1_12, spf_16, \
                         spf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = pb_y[k] * ssg0_10[k]
                  + f_0 * ssf_6[k]
                  - f_1 * pc_y[k] * ssg1_10[k];

        t_26[k] = f_2 * pc_z[k] * spf_16[k];

        t_27[k] = pb_y[k] * ssg0_12[k]
                  + f_3 * ssf_8[k]
                  - f_1 * pc_y[k] * ssg1_12[k];

        t_28[k] = f_4 * ssf_9[k]
                  + f_2 * pc_y[k] * spf_19[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pb_y, pb_z, pc_y, pc_z, ssg0_0, ssg0_14, \
                         ssf_0, ssg1_0, ssg1_14, spf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pb_y[k] * ssg0_14[k]
                  - f_1 * pc_y[k] * ssg1_14[k];

        t_30[k] = pb_z[k] * ssg0_0[k]
                  - f_1 * pc_z[k] * ssg1_0[k];

        t_31[k] = f_2 * pc_y[k] * spf_20[k];

        t_32[k] = f_4 * ssf_0[k]
                  + f_2 * pc_z[k] * spf_20[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pb_z, pc_x, pc_y, pc_z, ssg0_3, ssg1_3, \
                         spd0_17, spd1_17, spf_22, spf_25, spf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pb_z[k] * ssg0_3[k]
                  - f_1 * pc_z[k] * ssg1_3[k];

        t_34[k] = f_2 * pc_y[k] * spf_22[k];

        t_35[k] = f_5 * spd0_17[k]
                  - f_6 * spd1_17[k]
                  + f_2 * pc_x[k] * spf_25[k];

        t_36[k] = f_2 * pc_x[k] * spf_26[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, t_41, pb_z, pc_x, pc_z, ssg0_10, ssf_6, \
                         ssg1_10, spf_26, spf_27, spf_28, spf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_2 * pc_x[k] * spf_27[k];

        t_38[k] = f_2 * pc_x[k] * spf_28[k];

        t_39[k] = f_2 * pc_x[k] * spf_29[k];

        t_40[k] = pb_z[k] * ssg0_10[k]
                  - f_1 * pc_z[k] * ssg1_10[k];

        t_41[k] = f_4 * ssf_6[k]
                  + f_2 * pc_z[k] * spf_26[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pb_z, pc_y, pc_z, ssg0_14, ssf_9, ssg1_14, spd0_17, \
                         spd1_17, spf_28, spf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_5 * spd0_17[k]
                  - f_6 * spd1_17[k]
                  + f_2 * pc_y[k] * spf_28[k];

        t_43[k] = f_2 * pc_y[k] * spf_29[k];

        t_44[k] = pb_z[k] * ssg0_14[k]
                  + f_0 * ssf_9[k]
                  - f_1 * pc_z[k] * ssg1_14[k];
    }
}

}  // namespace simdt3ceri
