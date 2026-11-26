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


#include "SimdKineticEnergyVrrRecSI.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_prim_si_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pb,
                                 const size_t sg_s, const size_t si_s, const size_t sg,
                                 const size_t sh, const size_t ncols, const double alpha,
                                 const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.0 * alpha / p;
    const auto f_1 = 2.0 * alpha * beta / p;
    const auto f_2 = 2.5 / p;
    const auto f_3 = 3.0 * alpha / p;
    const auto f_4 = 1.5 / p;
    const auto f_5 = 2.0 * alpha / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = alpha / p;
    const auto f_8 = 0.5 / p;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *sg_s_0 = buffer.data(sg_s + 0);
    const auto *sg_s_3 = buffer.data(sg_s + 3);
    const auto *sg_s_5 = buffer.data(sg_s + 5);
    const auto *sg_s_6 = buffer.data(sg_s + 6);
    const auto *sg_s_9 = buffer.data(sg_s + 9);
    const auto *sg_s_10 = buffer.data(sg_s + 10);
    const auto *sg_s_12 = buffer.data(sg_s + 12);
    const auto *sg_s_13 = buffer.data(sg_s + 13);
    const auto *sg_s_14 = buffer.data(sg_s + 14);

    const auto *si_s_0 = buffer.data(si_s + 0);
    const auto *si_s_1 = buffer.data(si_s + 1);
    const auto *si_s_2 = buffer.data(si_s + 2);
    const auto *si_s_3 = buffer.data(si_s + 3);
    const auto *si_s_4 = buffer.data(si_s + 4);
    const auto *si_s_5 = buffer.data(si_s + 5);
    const auto *si_s_6 = buffer.data(si_s + 6);
    const auto *si_s_7 = buffer.data(si_s + 7);
    const auto *si_s_8 = buffer.data(si_s + 8);
    const auto *si_s_9 = buffer.data(si_s + 9);
    const auto *si_s_10 = buffer.data(si_s + 10);
    const auto *si_s_11 = buffer.data(si_s + 11);
    const auto *si_s_12 = buffer.data(si_s + 12);
    const auto *si_s_13 = buffer.data(si_s + 13);
    const auto *si_s_14 = buffer.data(si_s + 14);
    const auto *si_s_15 = buffer.data(si_s + 15);
    const auto *si_s_16 = buffer.data(si_s + 16);
    const auto *si_s_17 = buffer.data(si_s + 17);
    const auto *si_s_18 = buffer.data(si_s + 18);
    const auto *si_s_19 = buffer.data(si_s + 19);
    const auto *si_s_20 = buffer.data(si_s + 20);
    const auto *si_s_21 = buffer.data(si_s + 21);
    const auto *si_s_22 = buffer.data(si_s + 22);
    const auto *si_s_23 = buffer.data(si_s + 23);
    const auto *si_s_24 = buffer.data(si_s + 24);
    const auto *si_s_25 = buffer.data(si_s + 25);
    const auto *si_s_26 = buffer.data(si_s + 26);
    const auto *si_s_27 = buffer.data(si_s + 27);

    const auto *sg_0 = buffer.data(sg + 0);
    const auto *sg_3 = buffer.data(sg + 3);
    const auto *sg_5 = buffer.data(sg + 5);
    const auto *sg_6 = buffer.data(sg + 6);
    const auto *sg_9 = buffer.data(sg + 9);
    const auto *sg_10 = buffer.data(sg + 10);
    const auto *sg_12 = buffer.data(sg + 12);
    const auto *sg_13 = buffer.data(sg + 13);
    const auto *sg_14 = buffer.data(sg + 14);

    const auto *sh_0 = buffer.data(sh + 0);
    const auto *sh_2 = buffer.data(sh + 2);
    const auto *sh_3 = buffer.data(sh + 3);
    const auto *sh_5 = buffer.data(sh + 5);
    const auto *sh_6 = buffer.data(sh + 6);
    const auto *sh_9 = buffer.data(sh + 9);
    const auto *sh_10 = buffer.data(sh + 10);
    const auto *sh_12 = buffer.data(sh + 12);
    const auto *sh_14 = buffer.data(sh + 14);
    const auto *sh_15 = buffer.data(sh + 15);
    const auto *sh_16 = buffer.data(sh + 16);
    const auto *sh_17 = buffer.data(sh + 17);
    const auto *sh_18 = buffer.data(sh + 18);
    const auto *sh_19 = buffer.data(sh + 19);
    const auto *sh_20 = buffer.data(sh + 20);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, sg_s_0, si_s_0, si_s_1, si_s_2, \
                         sg_0, sh_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -f_0 * sg_s_0[k]
                 + f_1 * si_s_0[k]
                 + f_2 * sg_0[k]
                 + pb_x[k] * sh_0[k];

        t_1[k] = f_1 * si_s_1[k]
                 + pb_y[k] * sh_0[k];

        t_2[k] = f_1 * si_s_2[k]
                 + pb_z[k] * sh_0[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, pb_y, sg_s_3, sg_s_5, si_s_3, si_s_4, si_s_5, \
                         sg_3, sg_5, sh_2, sh_3, sh_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = -f_3 * sg_s_3[k]
                 + f_1 * si_s_3[k]
                 + f_4 * sg_3[k]
                 + pb_x[k] * sh_3[k];

        t_4[k] = f_1 * si_s_4[k]
                 + pb_y[k] * sh_2[k];

        t_5[k] = -f_3 * sg_s_5[k]
                 + f_1 * si_s_5[k]
                 + f_4 * sg_5[k]
                 + pb_x[k] * sh_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, pb_x, pb_y, pb_z, sg_s_6, si_s_6, si_s_7, si_s_8, \
                         sg_6, sh_3, sh_5, sh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_5 * sg_s_6[k]
                 + f_1 * si_s_6[k]
                 + f_6 * sg_6[k]
                 + pb_x[k] * sh_6[k];

        t_7[k] = f_1 * si_s_7[k]
                 + pb_z[k] * sh_3[k];

        t_8[k] = f_1 * si_s_8[k]
                 + pb_y[k] * sh_5[k];
    }

#pragma omp simd aligned(t_9, t_10, t_11, pb_x, pb_z, sg_s_9, sg_s_10, si_s_9, si_s_10, \
                         si_s_11, sg_9, sg_10, sh_6, sh_9, sh_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = -f_5 * sg_s_9[k]
                 + f_1 * si_s_9[k]
                 + f_6 * sg_9[k]
                 + pb_x[k] * sh_9[k];

        t_10[k] = -f_7 * sg_s_10[k]
                  + f_1 * si_s_10[k]
                  + f_8 * sg_10[k]
                  + pb_x[k] * sh_10[k];

        t_11[k] = f_1 * si_s_11[k]
                  + pb_z[k] * sh_6[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, pb_x, pb_y, sg_s_12, sg_s_14, si_s_12, si_s_13, \
                         si_s_14, sg_12, sg_14, sh_9, sh_12, sh_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -f_7 * sg_s_12[k]
                  + f_1 * si_s_12[k]
                  + f_8 * sg_12[k]
                  + pb_x[k] * sh_12[k];

        t_13[k] = f_1 * si_s_13[k]
                  + pb_y[k] * sh_9[k];

        t_14[k] = -f_7 * sg_s_14[k]
                  + f_1 * si_s_14[k]
                  + f_8 * sg_14[k]
                  + pb_x[k] * sh_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, pb_x, si_s_15, si_s_16, si_s_17, \
                         si_s_18, si_s_19, sh_15, sh_16, sh_17, sh_18, \
                         sh_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_1 * si_s_15[k]
                  + pb_x[k] * sh_15[k];

        t_16[k] = f_1 * si_s_16[k]
                  + pb_x[k] * sh_16[k];

        t_17[k] = f_1 * si_s_17[k]
                  + pb_x[k] * sh_17[k];

        t_18[k] = f_1 * si_s_18[k]
                  + pb_x[k] * sh_18[k];

        t_19[k] = f_1 * si_s_19[k]
                  + pb_x[k] * sh_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pb_x, pb_y, pb_z, sg_s_10, si_s_20, si_s_21, \
                         si_s_22, sg_10, sh_15, sh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_1 * si_s_20[k]
                  + pb_x[k] * sh_20[k];

        t_21[k] = -f_0 * sg_s_10[k]
                  + f_1 * si_s_21[k]
                  + f_2 * sg_10[k]
                  + pb_y[k] * sh_15[k];

        t_22[k] = f_1 * si_s_22[k]
                  + pb_z[k] * sh_15[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pb_y, sg_s_12, sg_s_13, sg_s_14, si_s_23, si_s_24, \
                         si_s_25, sg_12, sg_13, sg_14, sh_17, sh_18, \
                         sh_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = -f_3 * sg_s_12[k]
                  + f_1 * si_s_23[k]
                  + f_4 * sg_12[k]
                  + pb_y[k] * sh_17[k];

        t_24[k] = -f_5 * sg_s_13[k]
                  + f_1 * si_s_24[k]
                  + f_6 * sg_13[k]
                  + pb_y[k] * sh_18[k];

        t_25[k] = -f_7 * sg_s_14[k]
                  + f_1 * si_s_25[k]
                  + f_8 * sg_14[k]
                  + pb_y[k] * sh_19[k];
    }

#pragma omp simd aligned(t_26, t_27, pb_y, pb_z, sg_s_14, si_s_26, si_s_27, sg_14, \
                         sh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_1 * si_s_26[k]
                  + pb_y[k] * sh_20[k];

        t_27[k] = -f_0 * sg_s_14[k]
                  + f_1 * si_s_27[k]
                  + f_2 * sg_14[k]
                  + pb_z[k] * sh_20[k];
    }
}

}  // namespace simdkin
