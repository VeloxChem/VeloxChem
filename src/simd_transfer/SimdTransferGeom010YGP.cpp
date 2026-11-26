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


#include "SimdTransferGeom010YGP.hpp"

#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
compute_hrr_geom_010y_gp_out_of_second(CSimdMatrix &buffer, const CSimdMatrix &coordinates,
                                       const size_t target, const size_t fp_1, const size_t fp_0,
                                       const size_t fd_1, const size_t ncomps,
                                       const size_t nmax) -> void
{
    // NOTE: what the block carries below the pair reaches this step as a count.
    // The coefficients are powers of the separation and name no index of it, so
    // it is stepped over rather than known.

    for (size_t c = 0; c < ncomps; c++)
    {
        auto *t_0 = buffer.data(target + 0 * ncomps + c);
        auto *t_1 = buffer.data(target + 1 * ncomps + c);
        auto *t_2 = buffer.data(target + 2 * ncomps + c);
        auto *t_3 = buffer.data(target + 3 * ncomps + c);
        auto *t_4 = buffer.data(target + 4 * ncomps + c);
        auto *t_5 = buffer.data(target + 5 * ncomps + c);
        auto *t_6 = buffer.data(target + 6 * ncomps + c);
        auto *t_7 = buffer.data(target + 7 * ncomps + c);
        auto *t_8 = buffer.data(target + 8 * ncomps + c);
        auto *t_9 = buffer.data(target + 9 * ncomps + c);
        auto *t_10 = buffer.data(target + 10 * ncomps + c);
        auto *t_11 = buffer.data(target + 11 * ncomps + c);
        auto *t_12 = buffer.data(target + 12 * ncomps + c);
        auto *t_13 = buffer.data(target + 13 * ncomps + c);
        auto *t_14 = buffer.data(target + 14 * ncomps + c);
        auto *t_15 = buffer.data(target + 15 * ncomps + c);
        auto *t_16 = buffer.data(target + 16 * ncomps + c);
        auto *t_17 = buffer.data(target + 17 * ncomps + c);
        auto *t_18 = buffer.data(target + 18 * ncomps + c);
        auto *t_19 = buffer.data(target + 19 * ncomps + c);
        auto *t_20 = buffer.data(target + 20 * ncomps + c);
        auto *t_21 = buffer.data(target + 21 * ncomps + c);
        auto *t_22 = buffer.data(target + 22 * ncomps + c);
        auto *t_23 = buffer.data(target + 23 * ncomps + c);
        auto *t_24 = buffer.data(target + 24 * ncomps + c);
        auto *t_25 = buffer.data(target + 25 * ncomps + c);
        auto *t_26 = buffer.data(target + 26 * ncomps + c);
        auto *t_27 = buffer.data(target + 27 * ncomps + c);
        auto *t_28 = buffer.data(target + 28 * ncomps + c);
        auto *t_29 = buffer.data(target + 29 * ncomps + c);
        auto *t_30 = buffer.data(target + 30 * ncomps + c);
        auto *t_31 = buffer.data(target + 31 * ncomps + c);
        auto *t_32 = buffer.data(target + 32 * ncomps + c);
        auto *t_33 = buffer.data(target + 33 * ncomps + c);
        auto *t_34 = buffer.data(target + 34 * ncomps + c);
        auto *t_35 = buffer.data(target + 35 * ncomps + c);
        auto *t_36 = buffer.data(target + 36 * ncomps + c);
        auto *t_37 = buffer.data(target + 37 * ncomps + c);
        auto *t_38 = buffer.data(target + 38 * ncomps + c);
        auto *t_39 = buffer.data(target + 39 * ncomps + c);
        auto *t_40 = buffer.data(target + 40 * ncomps + c);
        auto *t_41 = buffer.data(target + 41 * ncomps + c);
        auto *t_42 = buffer.data(target + 42 * ncomps + c);
        auto *t_43 = buffer.data(target + 43 * ncomps + c);
        auto *t_44 = buffer.data(target + 44 * ncomps + c);

        const auto *ab_x = coordinates.data(6);
        const auto *ab_y = coordinates.data(7);
        const auto *ab_z = coordinates.data(8);

        const auto *fp_1_0 = buffer.data(fp_1 + 0 * ncomps + c);
        const auto *fp_1_1 = buffer.data(fp_1 + 1 * ncomps + c);
        const auto *fp_1_2 = buffer.data(fp_1 + 2 * ncomps + c);
        const auto *fp_1_3 = buffer.data(fp_1 + 3 * ncomps + c);
        const auto *fp_1_4 = buffer.data(fp_1 + 4 * ncomps + c);
        const auto *fp_1_5 = buffer.data(fp_1 + 5 * ncomps + c);
        const auto *fp_1_6 = buffer.data(fp_1 + 6 * ncomps + c);
        const auto *fp_1_7 = buffer.data(fp_1 + 7 * ncomps + c);
        const auto *fp_1_8 = buffer.data(fp_1 + 8 * ncomps + c);
        const auto *fp_1_9 = buffer.data(fp_1 + 9 * ncomps + c);
        const auto *fp_1_10 = buffer.data(fp_1 + 10 * ncomps + c);
        const auto *fp_1_11 = buffer.data(fp_1 + 11 * ncomps + c);
        const auto *fp_1_12 = buffer.data(fp_1 + 12 * ncomps + c);
        const auto *fp_1_13 = buffer.data(fp_1 + 13 * ncomps + c);
        const auto *fp_1_14 = buffer.data(fp_1 + 14 * ncomps + c);
        const auto *fp_1_15 = buffer.data(fp_1 + 15 * ncomps + c);
        const auto *fp_1_16 = buffer.data(fp_1 + 16 * ncomps + c);
        const auto *fp_1_17 = buffer.data(fp_1 + 17 * ncomps + c);
        const auto *fp_1_18 = buffer.data(fp_1 + 18 * ncomps + c);
        const auto *fp_1_19 = buffer.data(fp_1 + 19 * ncomps + c);
        const auto *fp_1_20 = buffer.data(fp_1 + 20 * ncomps + c);
        const auto *fp_1_21 = buffer.data(fp_1 + 21 * ncomps + c);
        const auto *fp_1_22 = buffer.data(fp_1 + 22 * ncomps + c);
        const auto *fp_1_23 = buffer.data(fp_1 + 23 * ncomps + c);
        const auto *fp_1_24 = buffer.data(fp_1 + 24 * ncomps + c);
        const auto *fp_1_25 = buffer.data(fp_1 + 25 * ncomps + c);
        const auto *fp_1_26 = buffer.data(fp_1 + 26 * ncomps + c);
        const auto *fp_1_27 = buffer.data(fp_1 + 27 * ncomps + c);
        const auto *fp_1_28 = buffer.data(fp_1 + 28 * ncomps + c);
        const auto *fp_1_29 = buffer.data(fp_1 + 29 * ncomps + c);

        const auto *fp_0_18 = buffer.data(fp_0 + 18 * ncomps + c);
        const auto *fp_0_19 = buffer.data(fp_0 + 19 * ncomps + c);
        const auto *fp_0_20 = buffer.data(fp_0 + 20 * ncomps + c);
        const auto *fp_0_21 = buffer.data(fp_0 + 21 * ncomps + c);
        const auto *fp_0_22 = buffer.data(fp_0 + 22 * ncomps + c);
        const auto *fp_0_23 = buffer.data(fp_0 + 23 * ncomps + c);
        const auto *fp_0_24 = buffer.data(fp_0 + 24 * ncomps + c);
        const auto *fp_0_25 = buffer.data(fp_0 + 25 * ncomps + c);
        const auto *fp_0_26 = buffer.data(fp_0 + 26 * ncomps + c);
        const auto *fp_0_27 = buffer.data(fp_0 + 27 * ncomps + c);
        const auto *fp_0_28 = buffer.data(fp_0 + 28 * ncomps + c);
        const auto *fp_0_29 = buffer.data(fp_0 + 29 * ncomps + c);

        const auto *fd_1_0 = buffer.data(fd_1 + 0 * ncomps + c);
        const auto *fd_1_1 = buffer.data(fd_1 + 1 * ncomps + c);
        const auto *fd_1_2 = buffer.data(fd_1 + 2 * ncomps + c);
        const auto *fd_1_6 = buffer.data(fd_1 + 6 * ncomps + c);
        const auto *fd_1_7 = buffer.data(fd_1 + 7 * ncomps + c);
        const auto *fd_1_8 = buffer.data(fd_1 + 8 * ncomps + c);
        const auto *fd_1_12 = buffer.data(fd_1 + 12 * ncomps + c);
        const auto *fd_1_13 = buffer.data(fd_1 + 13 * ncomps + c);
        const auto *fd_1_14 = buffer.data(fd_1 + 14 * ncomps + c);
        const auto *fd_1_18 = buffer.data(fd_1 + 18 * ncomps + c);
        const auto *fd_1_19 = buffer.data(fd_1 + 19 * ncomps + c);
        const auto *fd_1_20 = buffer.data(fd_1 + 20 * ncomps + c);
        const auto *fd_1_24 = buffer.data(fd_1 + 24 * ncomps + c);
        const auto *fd_1_25 = buffer.data(fd_1 + 25 * ncomps + c);
        const auto *fd_1_26 = buffer.data(fd_1 + 26 * ncomps + c);
        const auto *fd_1_30 = buffer.data(fd_1 + 30 * ncomps + c);
        const auto *fd_1_31 = buffer.data(fd_1 + 31 * ncomps + c);
        const auto *fd_1_32 = buffer.data(fd_1 + 32 * ncomps + c);
        const auto *fd_1_36 = buffer.data(fd_1 + 36 * ncomps + c);
        const auto *fd_1_37 = buffer.data(fd_1 + 37 * ncomps + c);
        const auto *fd_1_38 = buffer.data(fd_1 + 38 * ncomps + c);
        const auto *fd_1_39 = buffer.data(fd_1 + 39 * ncomps + c);
        const auto *fd_1_40 = buffer.data(fd_1 + 40 * ncomps + c);
        const auto *fd_1_42 = buffer.data(fd_1 + 42 * ncomps + c);
        const auto *fd_1_43 = buffer.data(fd_1 + 43 * ncomps + c);
        const auto *fd_1_44 = buffer.data(fd_1 + 44 * ncomps + c);
        const auto *fd_1_45 = buffer.data(fd_1 + 45 * ncomps + c);
        const auto *fd_1_46 = buffer.data(fd_1 + 46 * ncomps + c);
        const auto *fd_1_48 = buffer.data(fd_1 + 48 * ncomps + c);
        const auto *fd_1_49 = buffer.data(fd_1 + 49 * ncomps + c);
        const auto *fd_1_50 = buffer.data(fd_1 + 50 * ncomps + c);
        const auto *fd_1_51 = buffer.data(fd_1 + 51 * ncomps + c);
        const auto *fd_1_52 = buffer.data(fd_1 + 52 * ncomps + c);
        const auto *fd_1_54 = buffer.data(fd_1 + 54 * ncomps + c);
        const auto *fd_1_55 = buffer.data(fd_1 + 55 * ncomps + c);
        const auto *fd_1_56 = buffer.data(fd_1 + 56 * ncomps + c);
        const auto *fd_1_57 = buffer.data(fd_1 + 57 * ncomps + c);
        const auto *fd_1_58 = buffer.data(fd_1 + 58 * ncomps + c);
        const auto *fd_1_59 = buffer.data(fd_1 + 59 * ncomps + c);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, ab_x, fp_1_0, fp_1_1, fp_1_2, fp_1_3, \
                         fp_1_4, fd_1_0, fd_1_1, fd_1_2, fd_1_6, \
                         fd_1_7 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_0[k] = -ab_x[k] * fp_1_0[k]
                     + fd_1_0[k];

            t_1[k] = -ab_x[k] * fp_1_1[k]
                     + fd_1_1[k];

            t_2[k] = -ab_x[k] * fp_1_2[k]
                     + fd_1_2[k];

            t_3[k] = -ab_x[k] * fp_1_3[k]
                     + fd_1_6[k];

            t_4[k] = -ab_x[k] * fp_1_4[k]
                     + fd_1_7[k];
        }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, ab_x, fp_1_5, fp_1_6, fp_1_7, fp_1_8, \
                         fp_1_9, fd_1_8, fd_1_12, fd_1_13, fd_1_14, \
                         fd_1_18 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_5[k] = -ab_x[k] * fp_1_5[k]
                     + fd_1_8[k];

            t_6[k] = -ab_x[k] * fp_1_6[k]
                     + fd_1_12[k];

            t_7[k] = -ab_x[k] * fp_1_7[k]
                     + fd_1_13[k];

            t_8[k] = -ab_x[k] * fp_1_8[k]
                     + fd_1_14[k];

            t_9[k] = -ab_x[k] * fp_1_9[k]
                     + fd_1_18[k];
        }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, ab_x, fp_1_10, fp_1_11, fp_1_12, \
                         fp_1_13, fp_1_14, fd_1_19, fd_1_20, fd_1_24, fd_1_25, \
                         fd_1_26 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_10[k] = -ab_x[k] * fp_1_10[k]
                      + fd_1_19[k];

            t_11[k] = -ab_x[k] * fp_1_11[k]
                      + fd_1_20[k];

            t_12[k] = -ab_x[k] * fp_1_12[k]
                      + fd_1_24[k];

            t_13[k] = -ab_x[k] * fp_1_13[k]
                      + fd_1_25[k];

            t_14[k] = -ab_x[k] * fp_1_14[k]
                      + fd_1_26[k];
        }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, ab_x, fp_1_15, fp_1_16, fp_1_17, \
                         fp_1_18, fp_1_19, fd_1_30, fd_1_31, fd_1_32, fd_1_36, \
                         fd_1_37 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_15[k] = -ab_x[k] * fp_1_15[k]
                      + fd_1_30[k];

            t_16[k] = -ab_x[k] * fp_1_16[k]
                      + fd_1_31[k];

            t_17[k] = -ab_x[k] * fp_1_17[k]
                      + fd_1_32[k];

            t_18[k] = -ab_x[k] * fp_1_18[k]
                      + fd_1_36[k];

            t_19[k] = -ab_x[k] * fp_1_19[k]
                      + fd_1_37[k];
        }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, ab_x, fp_1_20, fp_1_21, fp_1_22, \
                         fp_1_23, fp_1_24, fd_1_38, fd_1_42, fd_1_43, fd_1_44, \
                         fd_1_48 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_20[k] = -ab_x[k] * fp_1_20[k]
                      + fd_1_38[k];

            t_21[k] = -ab_x[k] * fp_1_21[k]
                      + fd_1_42[k];

            t_22[k] = -ab_x[k] * fp_1_22[k]
                      + fd_1_43[k];

            t_23[k] = -ab_x[k] * fp_1_23[k]
                      + fd_1_44[k];

            t_24[k] = -ab_x[k] * fp_1_24[k]
                      + fd_1_48[k];
        }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, ab_x, fp_1_25, fp_1_26, fp_1_27, \
                         fp_1_28, fp_1_29, fd_1_49, fd_1_50, fd_1_54, fd_1_55, \
                         fd_1_56 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_25[k] = -ab_x[k] * fp_1_25[k]
                      + fd_1_49[k];

            t_26[k] = -ab_x[k] * fp_1_26[k]
                      + fd_1_50[k];

            t_27[k] = -ab_x[k] * fp_1_27[k]
                      + fd_1_54[k];

            t_28[k] = -ab_x[k] * fp_1_28[k]
                      + fd_1_55[k];

            t_29[k] = -ab_x[k] * fp_1_29[k]
                      + fd_1_56[k];
        }

#pragma omp simd aligned(t_30, t_31, t_32, ab_y, fp_1_18, fp_1_19, fp_1_20, fp_0_18, fp_0_19, \
                         fp_0_20, fd_1_37, fd_1_39, fd_1_40 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_30[k] = -ab_y[k] * fp_1_18[k]
                      + fp_0_18[k]
                      + fd_1_37[k];

            t_31[k] = -ab_y[k] * fp_1_19[k]
                      + fp_0_19[k]
                      + fd_1_39[k];

            t_32[k] = -ab_y[k] * fp_1_20[k]
                      + fp_0_20[k]
                      + fd_1_40[k];
        }

#pragma omp simd aligned(t_33, t_34, t_35, ab_y, fp_1_21, fp_1_22, fp_1_23, fp_0_21, fp_0_22, \
                         fp_0_23, fd_1_43, fd_1_45, fd_1_46 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_33[k] = -ab_y[k] * fp_1_21[k]
                      + fp_0_21[k]
                      + fd_1_43[k];

            t_34[k] = -ab_y[k] * fp_1_22[k]
                      + fp_0_22[k]
                      + fd_1_45[k];

            t_35[k] = -ab_y[k] * fp_1_23[k]
                      + fp_0_23[k]
                      + fd_1_46[k];
        }

#pragma omp simd aligned(t_36, t_37, t_38, ab_y, fp_1_24, fp_1_25, fp_1_26, fp_0_24, fp_0_25, \
                         fp_0_26, fd_1_49, fd_1_51, fd_1_52 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_36[k] = -ab_y[k] * fp_1_24[k]
                      + fp_0_24[k]
                      + fd_1_49[k];

            t_37[k] = -ab_y[k] * fp_1_25[k]
                      + fp_0_25[k]
                      + fd_1_51[k];

            t_38[k] = -ab_y[k] * fp_1_26[k]
                      + fp_0_26[k]
                      + fd_1_52[k];
        }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, ab_y, ab_z, fp_1_27, fp_1_28, fp_1_29, \
                         fp_0_27, fp_0_28, fp_0_29, fd_1_55, fd_1_56, fd_1_57, \
                         fd_1_58 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_39[k] = -ab_y[k] * fp_1_27[k]
                      + fp_0_27[k]
                      + fd_1_55[k];

            t_40[k] = -ab_y[k] * fp_1_28[k]
                      + fp_0_28[k]
                      + fd_1_57[k];

            t_41[k] = -ab_y[k] * fp_1_29[k]
                      + fp_0_29[k]
                      + fd_1_58[k];

            t_42[k] = -ab_z[k] * fp_1_27[k]
                      + fd_1_56[k];
        }

#pragma omp simd aligned(t_43, t_44, ab_z, fp_1_28, fp_1_29, fd_1_58, \
                         fd_1_59 : simd::cache_line_size())
        for (size_t k = 0; k < nmax; k++)
        {
            t_43[k] = -ab_z[k] * fp_1_28[k]
                      + fd_1_58[k];

            t_44[k] = -ab_z[k] * fp_1_29[k]
                      + fd_1_59[k];
        }
    }
}

}  // namespace simdtrf
