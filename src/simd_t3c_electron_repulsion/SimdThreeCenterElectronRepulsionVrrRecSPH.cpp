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


#include "SimdThreeCenterElectronRepulsionVrrRecSPH.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_sph_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pb, const size_t pc,
                                                   const size_t ssh0, const size_t ssg,
                                                   const size_t ssh1, const size_t spf0,
                                                   const size_t spf1, const size_t spg,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / q;
    const auto f_1 = gamma / q;
    const auto f_2 = p / q;
    const auto f_3 = 1.5 / q;
    const auto f_4 = 1.0 / q;
    const auto f_5 = 0.5 / q;
    const auto f_6 = 1.0 / gamma;
    const auto f_7 = p / (gamma * q);
    const auto f_8 = 0.5 / gamma;
    const auto f_9 = 0.5 * p / (gamma * q);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ssh0_0 = buffer.data(ssh0 + 0);
    const auto *ssh0_3 = buffer.data(ssh0 + 3);
    const auto *ssh0_5 = buffer.data(ssh0 + 5);
    const auto *ssh0_6 = buffer.data(ssh0 + 6);
    const auto *ssh0_9 = buffer.data(ssh0 + 9);
    const auto *ssh0_15 = buffer.data(ssh0 + 15);
    const auto *ssh0_17 = buffer.data(ssh0 + 17);
    const auto *ssh0_18 = buffer.data(ssh0 + 18);
    const auto *ssh0_20 = buffer.data(ssh0 + 20);

    const auto *ssg_0 = buffer.data(ssg + 0);
    const auto *ssg_2 = buffer.data(ssg + 2);
    const auto *ssg_3 = buffer.data(ssg + 3);
    const auto *ssg_5 = buffer.data(ssg + 5);
    const auto *ssg_6 = buffer.data(ssg + 6);
    const auto *ssg_9 = buffer.data(ssg + 9);
    const auto *ssg_10 = buffer.data(ssg + 10);
    const auto *ssg_11 = buffer.data(ssg + 11);
    const auto *ssg_12 = buffer.data(ssg + 12);
    const auto *ssg_13 = buffer.data(ssg + 13);
    const auto *ssg_14 = buffer.data(ssg + 14);

    const auto *ssh1_0 = buffer.data(ssh1 + 0);
    const auto *ssh1_3 = buffer.data(ssh1 + 3);
    const auto *ssh1_5 = buffer.data(ssh1 + 5);
    const auto *ssh1_6 = buffer.data(ssh1 + 6);
    const auto *ssh1_9 = buffer.data(ssh1 + 9);
    const auto *ssh1_15 = buffer.data(ssh1 + 15);
    const auto *ssh1_17 = buffer.data(ssh1 + 17);
    const auto *ssh1_18 = buffer.data(ssh1 + 18);
    const auto *ssh1_20 = buffer.data(ssh1 + 20);

    const auto *spf0_13 = buffer.data(spf0 + 13);
    const auto *spf0_16 = buffer.data(spf0 + 16);
    const auto *spf0_25 = buffer.data(spf0 + 25);
    const auto *spf0_28 = buffer.data(spf0 + 28);
    const auto *spf0_29 = buffer.data(spf0 + 29);

    const auto *spf1_13 = buffer.data(spf1 + 13);
    const auto *spf1_16 = buffer.data(spf1 + 16);
    const auto *spf1_25 = buffer.data(spf1 + 25);
    const auto *spf1_28 = buffer.data(spf1 + 28);
    const auto *spf1_29 = buffer.data(spf1 + 29);

    const auto *spg_0 = buffer.data(spg + 0);
    const auto *spg_2 = buffer.data(spg + 2);
    const auto *spg_3 = buffer.data(spg + 3);
    const auto *spg_5 = buffer.data(spg + 5);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pc_x, pc_y, pc_z, ssh0_0, ssh0_3, ssg_0, \
                         ssg_3, ssh1_0, ssh1_3, spg_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = pb_x[k] * ssh0_0[k]
                 + f_0 * ssg_0[k]
                 - f_1 * pc_x[k] * ssh1_0[k];

        t_1[k] = f_2 * pc_y[k] * spg_0[k];

        t_2[k] = f_2 * pc_z[k] * spg_0[k];

        t_3[k] = pb_x[k] * ssh0_3[k]
                 + f_3 * ssg_3[k]
                 - f_1 * pc_x[k] * ssh1_3[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_x, pc_x, pc_y, pc_z, ssh0_5, ssh0_6, ssg_5, \
                         ssg_6, ssh1_5, ssh1_6, spg_2, spg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * pc_y[k] * spg_2[k];

        t_5[k] = pb_x[k] * ssh0_5[k]
                 + f_3 * ssg_5[k]
                 - f_1 * pc_x[k] * ssh1_5[k];

        t_6[k] = pb_x[k] * ssh0_6[k]
                 + f_4 * ssg_6[k]
                 - f_1 * pc_x[k] * ssh1_6[k];

        t_7[k] = f_2 * pc_z[k] * spg_3[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, t_11, pb_x, pc_x, pc_y, ssh0_9, ssg_9, ssg_10, \
                         ssg_11, ssh1_9, spg_5, spg_10, spg_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_2 * pc_y[k] * spg_5[k];

        t_9[k] = pb_x[k] * ssh0_9[k]
                 + f_4 * ssg_9[k]
                 - f_1 * pc_x[k] * ssh1_9[k];

        t_10[k] = f_5 * ssg_10[k]
                  + f_2 * pc_x[k] * spg_10[k];

        t_11[k] = f_5 * ssg_11[k]
                  + f_2 * pc_x[k] * spg_11[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, pb_x, pc_x, ssh0_15, ssg_12, ssg_13, ssg_14, \
                         ssh1_15, spg_12, spg_13, spg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = f_5 * ssg_12[k]
                  + f_2 * pc_x[k] * spg_12[k];

        t_13[k] = f_5 * ssg_13[k]
                  + f_2 * pc_x[k] * spg_13[k];

        t_14[k] = f_5 * ssg_14[k]
                  + f_2 * pc_x[k] * spg_14[k];

        t_15[k] = pb_x[k] * ssh0_15[k]
                  - f_1 * pc_x[k] * ssh1_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pb_x, pc_x, pc_y, pc_z, ssh0_17, ssh0_18, \
                         ssh1_17, ssh1_18, spg_10, spg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_2 * pc_z[k] * spg_10[k];

        t_17[k] = pb_x[k] * ssh0_17[k]
                  - f_1 * pc_x[k] * ssh1_17[k];

        t_18[k] = pb_x[k] * ssh0_18[k]
                  - f_1 * pc_x[k] * ssh1_18[k];

        t_19[k] = f_2 * pc_y[k] * spg_14[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pb_x, pb_y, pc_x, pc_y, pc_z, ssh0_0, \
                         ssh0_20, ssg_0, ssh1_0, ssh1_20, spg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pb_x[k] * ssh0_20[k]
                  - f_1 * pc_x[k] * ssh1_20[k];

        t_21[k] = pb_y[k] * ssh0_0[k]
                  - f_1 * pc_y[k] * ssh1_0[k];

        t_22[k] = f_5 * ssg_0[k]
                  + f_2 * pc_y[k] * spg_15[k];

        t_23[k] = f_2 * pc_z[k] * spg_15[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, pb_y, pc_x, pc_y, ssh0_5, ssg_2, ssh1_5, spf0_13, \
                         spf1_13, spg_17, spg_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_6 * spf0_13[k]
                  - f_7 * spf1_13[k]
                  + f_2 * pc_x[k] * spg_18[k];

        t_25[k] = f_5 * ssg_2[k]
                  + f_2 * pc_y[k] * spg_17[k];

        t_26[k] = pb_y[k] * ssh0_5[k]
                  - f_1 * pc_y[k] * ssh1_5[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, t_30, pb_y, pc_x, pc_y, pc_z, ssh0_9, ssg_5, \
                         ssh1_9, spf0_16, spf1_16, spg_18, spg_20, \
                         spg_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_8 * spf0_16[k]
                  - f_9 * spf1_16[k]
                  + f_2 * pc_x[k] * spg_21[k];

        t_28[k] = f_2 * pc_z[k] * spg_18[k];

        t_29[k] = f_5 * ssg_5[k]
                  + f_2 * pc_y[k] * spg_20[k];

        t_30[k] = pb_y[k] * ssh0_9[k]
                  - f_1 * pc_y[k] * ssh1_9[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pc_x, spg_25, spg_26, spg_27, spg_28, \
                         spg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_2 * pc_x[k] * spg_25[k];

        t_32[k] = f_2 * pc_x[k] * spg_26[k];

        t_33[k] = f_2 * pc_x[k] * spg_27[k];

        t_34[k] = f_2 * pc_x[k] * spg_28[k];

        t_35[k] = f_2 * pc_x[k] * spg_29[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pb_y, pc_y, pc_z, ssh0_15, ssh0_17, ssg_10, ssg_12, \
                         ssh1_15, ssh1_17, spg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = pb_y[k] * ssh0_15[k]
                  + f_0 * ssg_10[k]
                  - f_1 * pc_y[k] * ssh1_15[k];

        t_37[k] = f_2 * pc_z[k] * spg_25[k];

        t_38[k] = pb_y[k] * ssh0_17[k]
                  + f_3 * ssg_12[k]
                  - f_1 * pc_y[k] * ssh1_17[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_y, pc_y, ssh0_18, ssh0_20, ssg_13, ssg_14, \
                         ssh1_18, ssh1_20, spg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = pb_y[k] * ssh0_18[k]
                  + f_4 * ssg_13[k]
                  - f_1 * pc_y[k] * ssh1_18[k];

        t_40[k] = f_5 * ssg_14[k]
                  + f_2 * pc_y[k] * spg_29[k];

        t_41[k] = pb_y[k] * ssh0_20[k]
                  - f_1 * pc_y[k] * ssh1_20[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, t_46, pb_z, pc_y, pc_z, ssh0_0, ssh0_3, \
                         ssg_0, ssh1_0, ssh1_3, spg_30, spg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pb_z[k] * ssh0_0[k]
                  - f_1 * pc_z[k] * ssh1_0[k];

        t_43[k] = f_2 * pc_y[k] * spg_30[k];

        t_44[k] = f_5 * ssg_0[k]
                  + f_2 * pc_z[k] * spg_30[k];

        t_45[k] = pb_z[k] * ssh0_3[k]
                  - f_1 * pc_z[k] * ssh1_3[k];

        t_46[k] = f_2 * pc_y[k] * spg_32[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pb_z, pc_x, pc_y, pc_z, ssh0_6, ssg_3, \
                         ssh1_6, spf0_25, spf1_25, spg_33, spg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_6 * spf0_25[k]
                  - f_7 * spf1_25[k]
                  + f_2 * pc_x[k] * spg_35[k];

        t_48[k] = pb_z[k] * ssh0_6[k]
                  - f_1 * pc_z[k] * ssh1_6[k];

        t_49[k] = f_5 * ssg_3[k]
                  + f_2 * pc_z[k] * spg_33[k];

        t_50[k] = f_2 * pc_y[k] * spg_35[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, t_56, pc_x, spf0_29, spf1_29, spg_39, \
                         spg_40, spg_41, spg_42, spg_43, spg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_8 * spf0_29[k]
                  - f_9 * spf1_29[k]
                  + f_2 * pc_x[k] * spg_39[k];

        t_52[k] = f_2 * pc_x[k] * spg_40[k];

        t_53[k] = f_2 * pc_x[k] * spg_41[k];

        t_54[k] = f_2 * pc_x[k] * spg_42[k];

        t_55[k] = f_2 * pc_x[k] * spg_43[k];

        t_56[k] = f_2 * pc_x[k] * spg_44[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, pb_z, pc_y, pc_z, ssh0_15, ssg_10, ssh1_15, \
                         spf0_28, spf1_28, spg_40, spg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = pb_z[k] * ssh0_15[k]
                  - f_1 * pc_z[k] * ssh1_15[k];

        t_58[k] = f_5 * ssg_10[k]
                  + f_2 * pc_z[k] * spg_40[k];

        t_59[k] = f_6 * spf0_28[k]
                  - f_7 * spf1_28[k]
                  + f_2 * pc_y[k] * spg_42[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, pb_z, pc_y, pc_z, ssh0_20, ssg_14, ssh1_20, \
                         spf0_29, spf1_29, spg_43, spg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_8 * spf0_29[k]
                  - f_9 * spf1_29[k]
                  + f_2 * pc_y[k] * spg_43[k];

        t_61[k] = f_2 * pc_y[k] * spg_44[k];

        t_62[k] = pb_z[k] * ssh0_20[k]
                  + f_0 * ssg_14[k]
                  - f_1 * pc_z[k] * ssh1_20[k];
    }
}

}  // namespace simdt3ceri
