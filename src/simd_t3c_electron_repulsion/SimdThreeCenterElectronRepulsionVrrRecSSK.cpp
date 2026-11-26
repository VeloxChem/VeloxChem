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


#include "SimdThreeCenterElectronRepulsionVrrRecSSK.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_ssk_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pc, const size_t ssh0,
                                                   const size_t ssh1, const size_t ssi,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / gamma;
    const auto f_1 = 3.0 * p / (gamma * q);
    const auto f_2 = p / q;
    const auto f_3 = 2.0 / gamma;
    const auto f_4 = 2.0 * p / (gamma * q);
    const auto f_5 = 1.5 / gamma;
    const auto f_6 = 1.5 * p / (gamma * q);
    const auto f_7 = 1.0 / gamma;
    const auto f_8 = p / (gamma * q);
    const auto f_9 = 0.5 / gamma;
    const auto f_10 = 0.5 * p / (gamma * q);

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ssh0_0 = buffer.data(ssh0 + 0);
    const auto *ssh0_3 = buffer.data(ssh0 + 3);
    const auto *ssh0_5 = buffer.data(ssh0 + 5);
    const auto *ssh0_6 = buffer.data(ssh0 + 6);
    const auto *ssh0_9 = buffer.data(ssh0 + 9);
    const auto *ssh0_10 = buffer.data(ssh0 + 10);
    const auto *ssh0_12 = buffer.data(ssh0 + 12);
    const auto *ssh0_14 = buffer.data(ssh0 + 14);
    const auto *ssh0_15 = buffer.data(ssh0 + 15);
    const auto *ssh0_17 = buffer.data(ssh0 + 17);
    const auto *ssh0_18 = buffer.data(ssh0 + 18);
    const auto *ssh0_19 = buffer.data(ssh0 + 19);
    const auto *ssh0_20 = buffer.data(ssh0 + 20);

    const auto *ssh1_0 = buffer.data(ssh1 + 0);
    const auto *ssh1_3 = buffer.data(ssh1 + 3);
    const auto *ssh1_5 = buffer.data(ssh1 + 5);
    const auto *ssh1_6 = buffer.data(ssh1 + 6);
    const auto *ssh1_9 = buffer.data(ssh1 + 9);
    const auto *ssh1_10 = buffer.data(ssh1 + 10);
    const auto *ssh1_12 = buffer.data(ssh1 + 12);
    const auto *ssh1_14 = buffer.data(ssh1 + 14);
    const auto *ssh1_15 = buffer.data(ssh1 + 15);
    const auto *ssh1_17 = buffer.data(ssh1 + 17);
    const auto *ssh1_18 = buffer.data(ssh1 + 18);
    const auto *ssh1_19 = buffer.data(ssh1 + 19);
    const auto *ssh1_20 = buffer.data(ssh1 + 20);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, pc_z, ssh0_0, ssh0_3, ssh1_0, \
                         ssh1_3, ssi_0, ssi_2, ssi_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ssh0_0[k]
                 - f_1 * ssh1_0[k]
                 + f_2 * pc_x[k] * ssi_0[k];

        t_1[k] = f_2 * pc_y[k] * ssi_0[k];

        t_2[k] = f_2 * pc_z[k] * ssi_0[k];

        t_3[k] = f_3 * ssh0_3[k]
                 - f_4 * ssh1_3[k]
                 + f_2 * pc_x[k] * ssi_3[k];

        t_4[k] = f_2 * pc_y[k] * ssi_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pc_x, pc_y, pc_z, ssh0_5, ssh0_6, ssh1_5, ssh1_6, \
                         ssi_3, ssi_5, ssi_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_3 * ssh0_5[k]
                 - f_4 * ssh1_5[k]
                 + f_2 * pc_x[k] * ssi_5[k];

        t_6[k] = f_5 * ssh0_6[k]
                 - f_6 * ssh1_6[k]
                 + f_2 * pc_x[k] * ssi_6[k];

        t_7[k] = f_2 * pc_z[k] * ssi_3[k];

        t_8[k] = f_2 * pc_y[k] * ssi_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pc_x, pc_z, ssh0_9, ssh0_10, ssh0_12, ssh1_9, \
                         ssh1_10, ssh1_12, ssi_6, ssi_9, ssi_10, \
                         ssi_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * ssh0_9[k]
                 - f_6 * ssh1_9[k]
                 + f_2 * pc_x[k] * ssi_9[k];

        t_10[k] = f_7 * ssh0_10[k]
                  - f_8 * ssh1_10[k]
                  + f_2 * pc_x[k] * ssi_10[k];

        t_11[k] = f_2 * pc_z[k] * ssi_6[k];

        t_12[k] = f_7 * ssh0_12[k]
                  - f_8 * ssh1_12[k]
                  + f_2 * pc_x[k] * ssi_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, pc_x, pc_y, pc_z, ssh0_14, ssh0_15, ssh1_14, \
                         ssh1_15, ssi_9, ssi_10, ssi_14, ssi_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_2 * pc_y[k] * ssi_9[k];

        t_14[k] = f_7 * ssh0_14[k]
                  - f_8 * ssh1_14[k]
                  + f_2 * pc_x[k] * ssi_14[k];

        t_15[k] = f_9 * ssh0_15[k]
                  - f_10 * ssh1_15[k]
                  + f_2 * pc_x[k] * ssi_15[k];

        t_16[k] = f_2 * pc_z[k] * ssi_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, t_20, pc_x, pc_y, ssh0_17, ssh0_18, ssh0_20, \
                         ssh1_17, ssh1_18, ssh1_20, ssi_14, ssi_17, ssi_18, \
                         ssi_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_9 * ssh0_17[k]
                  - f_10 * ssh1_17[k]
                  + f_2 * pc_x[k] * ssi_17[k];

        t_18[k] = f_9 * ssh0_18[k]
                  - f_10 * ssh1_18[k]
                  + f_2 * pc_x[k] * ssi_18[k];

        t_19[k] = f_2 * pc_y[k] * ssi_14[k];

        t_20[k] = f_9 * ssh0_20[k]
                  - f_10 * ssh1_20[k]
                  + f_2 * pc_x[k] * ssi_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, t_25, t_26, t_27, pc_x, ssi_21, ssi_22, \
                         ssi_23, ssi_24, ssi_25, ssi_26, ssi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_2 * pc_x[k] * ssi_21[k];

        t_22[k] = f_2 * pc_x[k] * ssi_22[k];

        t_23[k] = f_2 * pc_x[k] * ssi_23[k];

        t_24[k] = f_2 * pc_x[k] * ssi_24[k];

        t_25[k] = f_2 * pc_x[k] * ssi_25[k];

        t_26[k] = f_2 * pc_x[k] * ssi_26[k];

        t_27[k] = f_2 * pc_x[k] * ssi_27[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pc_y, pc_z, ssh0_15, ssh0_17, ssh0_18, \
                         ssh1_15, ssh1_17, ssh1_18, ssi_21, ssi_23, \
                         ssi_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = f_0 * ssh0_15[k]
                  - f_1 * ssh1_15[k]
                  + f_2 * pc_y[k] * ssi_21[k];

        t_29[k] = f_2 * pc_z[k] * ssi_21[k];

        t_30[k] = f_3 * ssh0_17[k]
                  - f_4 * ssh1_17[k]
                  + f_2 * pc_y[k] * ssi_23[k];

        t_31[k] = f_5 * ssh0_18[k]
                  - f_6 * ssh1_18[k]
                  + f_2 * pc_y[k] * ssi_24[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pc_y, pc_z, ssh0_19, ssh0_20, ssh1_19, \
                         ssh1_20, ssi_25, ssi_26, ssi_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_7 * ssh0_19[k]
                  - f_8 * ssh1_19[k]
                  + f_2 * pc_y[k] * ssi_25[k];

        t_33[k] = f_9 * ssh0_20[k]
                  - f_10 * ssh1_20[k]
                  + f_2 * pc_y[k] * ssi_26[k];

        t_34[k] = f_2 * pc_y[k] * ssi_27[k];

        t_35[k] = f_0 * ssh0_20[k]
                  - f_1 * ssh1_20[k]
                  + f_2 * pc_z[k] * ssi_27[k];
    }
}

}  // namespace simdt3ceri
