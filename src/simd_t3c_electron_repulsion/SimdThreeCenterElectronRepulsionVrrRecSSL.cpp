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


#include "SimdThreeCenterElectronRepulsionVrrRecSSL.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_ssl_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pc, const size_t ssi0,
                                                   const size_t ssi1, const size_t ssk,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.5 / gamma;
    const auto f_1 = 3.5 * p / (gamma * q);
    const auto f_2 = p / q;
    const auto f_3 = 2.5 / gamma;
    const auto f_4 = 2.5 * p / (gamma * q);
    const auto f_5 = 2.0 / gamma;
    const auto f_6 = 2.0 * p / (gamma * q);
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
    const auto *ssi0_15 = buffer.data(ssi0 + 15);
    const auto *ssi0_17 = buffer.data(ssi0 + 17);
    const auto *ssi0_18 = buffer.data(ssi0 + 18);
    const auto *ssi0_20 = buffer.data(ssi0 + 20);
    const auto *ssi0_21 = buffer.data(ssi0 + 21);
    const auto *ssi0_23 = buffer.data(ssi0 + 23);
    const auto *ssi0_24 = buffer.data(ssi0 + 24);
    const auto *ssi0_25 = buffer.data(ssi0 + 25);
    const auto *ssi0_26 = buffer.data(ssi0 + 26);
    const auto *ssi0_27 = buffer.data(ssi0 + 27);

    const auto *ssi1_0 = buffer.data(ssi1 + 0);
    const auto *ssi1_3 = buffer.data(ssi1 + 3);
    const auto *ssi1_5 = buffer.data(ssi1 + 5);
    const auto *ssi1_6 = buffer.data(ssi1 + 6);
    const auto *ssi1_9 = buffer.data(ssi1 + 9);
    const auto *ssi1_10 = buffer.data(ssi1 + 10);
    const auto *ssi1_12 = buffer.data(ssi1 + 12);
    const auto *ssi1_14 = buffer.data(ssi1 + 14);
    const auto *ssi1_15 = buffer.data(ssi1 + 15);
    const auto *ssi1_17 = buffer.data(ssi1 + 17);
    const auto *ssi1_18 = buffer.data(ssi1 + 18);
    const auto *ssi1_20 = buffer.data(ssi1 + 20);
    const auto *ssi1_21 = buffer.data(ssi1 + 21);
    const auto *ssi1_23 = buffer.data(ssi1 + 23);
    const auto *ssi1_24 = buffer.data(ssi1 + 24);
    const auto *ssi1_25 = buffer.data(ssi1 + 25);
    const auto *ssi1_26 = buffer.data(ssi1 + 26);
    const auto *ssi1_27 = buffer.data(ssi1 + 27);

    const auto *ssk_0 = buffer.data(ssk + 0);
    const auto *ssk_2 = buffer.data(ssk + 2);
    const auto *ssk_3 = buffer.data(ssk + 3);
    const auto *ssk_5 = buffer.data(ssk + 5);
    const auto *ssk_6 = buffer.data(ssk + 6);
    const auto *ssk_9 = buffer.data(ssk + 9);
    const auto *ssk_10 = buffer.data(ssk + 10);
    const auto *ssk_12 = buffer.data(ssk + 12);
    const auto *ssk_14 = buffer.data(ssk + 14);
    const auto *ssk_15 = buffer.data(ssk + 15);
    const auto *ssk_17 = buffer.data(ssk + 17);
    const auto *ssk_18 = buffer.data(ssk + 18);
    const auto *ssk_20 = buffer.data(ssk + 20);
    const auto *ssk_21 = buffer.data(ssk + 21);
    const auto *ssk_23 = buffer.data(ssk + 23);
    const auto *ssk_24 = buffer.data(ssk + 24);
    const auto *ssk_25 = buffer.data(ssk + 25);
    const auto *ssk_27 = buffer.data(ssk + 27);
    const auto *ssk_28 = buffer.data(ssk + 28);
    const auto *ssk_29 = buffer.data(ssk + 29);
    const auto *ssk_30 = buffer.data(ssk + 30);
    const auto *ssk_31 = buffer.data(ssk + 31);
    const auto *ssk_32 = buffer.data(ssk + 32);
    const auto *ssk_33 = buffer.data(ssk + 33);
    const auto *ssk_34 = buffer.data(ssk + 34);
    const auto *ssk_35 = buffer.data(ssk + 35);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, pc_z, ssi0_0, ssi0_3, ssi1_0, \
                         ssi1_3, ssk_0, ssk_2, ssk_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ssi0_0[k]
                 - f_1 * ssi1_0[k]
                 + f_2 * pc_x[k] * ssk_0[k];

        t_1[k] = f_2 * pc_y[k] * ssk_0[k];

        t_2[k] = f_2 * pc_z[k] * ssk_0[k];

        t_3[k] = f_3 * ssi0_3[k]
                 - f_4 * ssi1_3[k]
                 + f_2 * pc_x[k] * ssk_3[k];

        t_4[k] = f_2 * pc_y[k] * ssk_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pc_x, pc_y, pc_z, ssi0_5, ssi0_6, ssi1_5, ssi1_6, \
                         ssk_3, ssk_5, ssk_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_3 * ssi0_5[k]
                 - f_4 * ssi1_5[k]
                 + f_2 * pc_x[k] * ssk_5[k];

        t_6[k] = f_5 * ssi0_6[k]
                 - f_6 * ssi1_6[k]
                 + f_2 * pc_x[k] * ssk_6[k];

        t_7[k] = f_2 * pc_z[k] * ssk_3[k];

        t_8[k] = f_2 * pc_y[k] * ssk_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pc_x, pc_z, ssi0_9, ssi0_10, ssi0_12, ssi1_9, \
                         ssi1_10, ssi1_12, ssk_6, ssk_9, ssk_10, \
                         ssk_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * ssi0_9[k]
                 - f_6 * ssi1_9[k]
                 + f_2 * pc_x[k] * ssk_9[k];

        t_10[k] = f_7 * ssi0_10[k]
                  - f_8 * ssi1_10[k]
                  + f_2 * pc_x[k] * ssk_10[k];

        t_11[k] = f_2 * pc_z[k] * ssk_6[k];

        t_12[k] = f_7 * ssi0_12[k]
                  - f_8 * ssi1_12[k]
                  + f_2 * pc_x[k] * ssk_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pc_x, pc_y, pc_z, ssi0_14, ssi0_15, ssi1_14, \
                         ssi1_15, ssk_9, ssk_10, ssk_14, ssk_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_2 * pc_y[k] * ssk_9[k];

        t_14[k] = f_7 * ssi0_14[k]
                  - f_8 * ssi1_14[k]
                  + f_2 * pc_x[k] * ssk_14[k];

        t_15[k] = f_9 * ssi0_15[k]
                  - f_10 * ssi1_15[k]
                  + f_2 * pc_x[k] * ssk_15[k];

        t_16[k] = f_2 * pc_z[k] * ssk_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pc_x, pc_y, ssi0_17, ssi0_18, ssi0_20, \
                         ssi1_17, ssi1_18, ssi1_20, ssk_14, ssk_17, ssk_18, \
                         ssk_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_9 * ssi0_17[k]
                  - f_10 * ssi1_17[k]
                  + f_2 * pc_x[k] * ssk_17[k];

        t_18[k] = f_9 * ssi0_18[k]
                  - f_10 * ssi1_18[k]
                  + f_2 * pc_x[k] * ssk_18[k];

        t_19[k] = f_2 * pc_y[k] * ssk_14[k];

        t_20[k] = f_9 * ssi0_20[k]
                  - f_10 * ssi1_20[k]
                  + f_2 * pc_x[k] * ssk_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pc_x, pc_z, ssi0_21, ssi0_23, ssi0_24, \
                         ssi1_21, ssi1_23, ssi1_24, ssk_15, ssk_21, ssk_23, \
                         ssk_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_11 * ssi0_21[k]
                  - f_12 * ssi1_21[k]
                  + f_2 * pc_x[k] * ssk_21[k];

        t_22[k] = f_2 * pc_z[k] * ssk_15[k];

        t_23[k] = f_11 * ssi0_23[k]
                  - f_12 * ssi1_23[k]
                  + f_2 * pc_x[k] * ssk_23[k];

        t_24[k] = f_11 * ssi0_24[k]
                  - f_12 * ssi1_24[k]
                  + f_2 * pc_x[k] * ssk_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pc_x, pc_y, ssi0_25, ssi0_27, ssi1_25, \
                         ssi1_27, ssk_20, ssk_25, ssk_27, ssk_28, \
                         ssk_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_11 * ssi0_25[k]
                  - f_12 * ssi1_25[k]
                  + f_2 * pc_x[k] * ssk_25[k];

        t_26[k] = f_2 * pc_y[k] * ssk_20[k];

        t_27[k] = f_11 * ssi0_27[k]
                  - f_12 * ssi1_27[k]
                  + f_2 * pc_x[k] * ssk_27[k];

        t_28[k] = f_2 * pc_x[k] * ssk_28[k];

        t_29[k] = f_2 * pc_x[k] * ssk_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, t_35, pc_x, ssk_30, ssk_31, ssk_32, \
                         ssk_33, ssk_34, ssk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_2 * pc_x[k] * ssk_30[k];

        t_31[k] = f_2 * pc_x[k] * ssk_31[k];

        t_32[k] = f_2 * pc_x[k] * ssk_32[k];

        t_33[k] = f_2 * pc_x[k] * ssk_33[k];

        t_34[k] = f_2 * pc_x[k] * ssk_34[k];

        t_35[k] = f_2 * pc_x[k] * ssk_35[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pc_y, pc_z, ssi0_21, ssi0_23, ssi0_24, \
                         ssi1_21, ssi1_23, ssi1_24, ssk_28, ssk_30, \
                         ssk_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_0 * ssi0_21[k]
                  - f_1 * ssi1_21[k]
                  + f_2 * pc_y[k] * ssk_28[k];

        t_37[k] = f_2 * pc_z[k] * ssk_28[k];

        t_38[k] = f_3 * ssi0_23[k]
                  - f_4 * ssi1_23[k]
                  + f_2 * pc_y[k] * ssk_30[k];

        t_39[k] = f_5 * ssi0_24[k]
                  - f_6 * ssi1_24[k]
                  + f_2 * pc_y[k] * ssk_31[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pc_y, ssi0_25, ssi0_26, ssi0_27, ssi1_25, \
                         ssi1_26, ssi1_27, ssk_32, ssk_33, ssk_34, \
                         ssk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_7 * ssi0_25[k]
                  - f_8 * ssi1_25[k]
                  + f_2 * pc_y[k] * ssk_32[k];

        t_41[k] = f_9 * ssi0_26[k]
                  - f_10 * ssi1_26[k]
                  + f_2 * pc_y[k] * ssk_33[k];

        t_42[k] = f_11 * ssi0_27[k]
                  - f_12 * ssi1_27[k]
                  + f_2 * pc_y[k] * ssk_34[k];

        t_43[k] = f_2 * pc_y[k] * ssk_35[k];
    }

#pragma omp simd aligned(t_44, pc_z, ssi0_27, ssi1_27, ssk_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_0 * ssi0_27[k]
                  - f_1 * ssi1_27[k]
                  + f_2 * pc_z[k] * ssk_35[k];
    }
}

}  // namespace simdt3ceri
