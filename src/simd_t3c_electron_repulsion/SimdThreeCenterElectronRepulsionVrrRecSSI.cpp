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


#include "SimdThreeCenterElectronRepulsionVrrRecSSI.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_ssi_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pc, const size_t ssg0,
                                                   const size_t ssg1, const size_t ssh,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / gamma;
    const auto f_1 = 2.5 * p / (gamma * q);
    const auto f_2 = p / q;
    const auto f_3 = 1.5 / gamma;
    const auto f_4 = 1.5 * p / (gamma * q);
    const auto f_5 = 1.0 / gamma;
    const auto f_6 = p / (gamma * q);
    const auto f_7 = 0.5 / gamma;
    const auto f_8 = 0.5 * p / (gamma * q);

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ssg0_0 = buffer.data(ssg0 + 0);
    const auto *ssg0_3 = buffer.data(ssg0 + 3);
    const auto *ssg0_5 = buffer.data(ssg0 + 5);
    const auto *ssg0_6 = buffer.data(ssg0 + 6);
    const auto *ssg0_9 = buffer.data(ssg0 + 9);
    const auto *ssg0_10 = buffer.data(ssg0 + 10);
    const auto *ssg0_12 = buffer.data(ssg0 + 12);
    const auto *ssg0_13 = buffer.data(ssg0 + 13);
    const auto *ssg0_14 = buffer.data(ssg0 + 14);

    const auto *ssg1_0 = buffer.data(ssg1 + 0);
    const auto *ssg1_3 = buffer.data(ssg1 + 3);
    const auto *ssg1_5 = buffer.data(ssg1 + 5);
    const auto *ssg1_6 = buffer.data(ssg1 + 6);
    const auto *ssg1_9 = buffer.data(ssg1 + 9);
    const auto *ssg1_10 = buffer.data(ssg1 + 10);
    const auto *ssg1_12 = buffer.data(ssg1 + 12);
    const auto *ssg1_13 = buffer.data(ssg1 + 13);
    const auto *ssg1_14 = buffer.data(ssg1 + 14);

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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, pc_x, pc_y, pc_z, ssg0_0, ssg0_3, ssg1_0, \
                         ssg1_3, ssh_0, ssh_2, ssh_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ssg0_0[k]
                 - f_1 * ssg1_0[k]
                 + f_2 * pc_x[k] * ssh_0[k];

        t_1[k] = f_2 * pc_y[k] * ssh_0[k];

        t_2[k] = f_2 * pc_z[k] * ssh_0[k];

        t_3[k] = f_3 * ssg0_3[k]
                 - f_4 * ssg1_3[k]
                 + f_2 * pc_x[k] * ssh_3[k];

        t_4[k] = f_2 * pc_y[k] * ssh_2[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, pc_x, pc_y, pc_z, ssg0_5, ssg0_6, ssg1_5, ssg1_6, \
                         ssh_3, ssh_5, ssh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_3 * ssg0_5[k]
                 - f_4 * ssg1_5[k]
                 + f_2 * pc_x[k] * ssh_5[k];

        t_6[k] = f_5 * ssg0_6[k]
                 - f_6 * ssg1_6[k]
                 + f_2 * pc_x[k] * ssh_6[k];

        t_7[k] = f_2 * pc_z[k] * ssh_3[k];

        t_8[k] = f_2 * pc_y[k] * ssh_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, t_12, pc_x, pc_z, ssg0_9, ssg0_10, ssg0_12, ssg1_9, \
                         ssg1_10, ssg1_12, ssh_6, ssh_9, ssh_10, \
                         ssh_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * ssg0_9[k]
                 - f_6 * ssg1_9[k]
                 + f_2 * pc_x[k] * ssh_9[k];

        t_10[k] = f_7 * ssg0_10[k]
                  - f_8 * ssg1_10[k]
                  + f_2 * pc_x[k] * ssh_10[k];

        t_11[k] = f_2 * pc_z[k] * ssh_6[k];

        t_12[k] = f_7 * ssg0_12[k]
                  - f_8 * ssg1_12[k]
                  + f_2 * pc_x[k] * ssh_12[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, t_18, pc_x, pc_y, ssg0_14, ssg1_14, \
                         ssh_9, ssh_14, ssh_15, ssh_16, ssh_17, \
                         ssh_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_2 * pc_y[k] * ssh_9[k];

        t_14[k] = f_7 * ssg0_14[k]
                  - f_8 * ssg1_14[k]
                  + f_2 * pc_x[k] * ssh_14[k];

        t_15[k] = f_2 * pc_x[k] * ssh_15[k];

        t_16[k] = f_2 * pc_x[k] * ssh_16[k];

        t_17[k] = f_2 * pc_x[k] * ssh_17[k];

        t_18[k] = f_2 * pc_x[k] * ssh_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, t_23, pc_x, pc_y, pc_z, ssg0_10, ssg0_12, \
                         ssg1_10, ssg1_12, ssh_15, ssh_17, ssh_19, \
                         ssh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_2 * pc_x[k] * ssh_19[k];

        t_20[k] = f_2 * pc_x[k] * ssh_20[k];

        t_21[k] = f_0 * ssg0_10[k]
                  - f_1 * ssg1_10[k]
                  + f_2 * pc_y[k] * ssh_15[k];

        t_22[k] = f_2 * pc_z[k] * ssh_15[k];

        t_23[k] = f_3 * ssg0_12[k]
                  - f_4 * ssg1_12[k]
                  + f_2 * pc_y[k] * ssh_17[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pc_y, pc_z, ssg0_13, ssg0_14, ssg1_13, \
                         ssg1_14, ssh_18, ssh_19, ssh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_5 * ssg0_13[k]
                  - f_6 * ssg1_13[k]
                  + f_2 * pc_y[k] * ssh_18[k];

        t_25[k] = f_7 * ssg0_14[k]
                  - f_8 * ssg1_14[k]
                  + f_2 * pc_y[k] * ssh_19[k];

        t_26[k] = f_2 * pc_y[k] * ssh_20[k];

        t_27[k] = f_0 * ssg0_14[k]
                  - f_1 * ssg1_14[k]
                  + f_2 * pc_z[k] * ssh_20[k];
    }
}

}  // namespace simdt3ceri
