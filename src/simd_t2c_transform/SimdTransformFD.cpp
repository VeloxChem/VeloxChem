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


#include "SimdTransformFD.hpp"

#include <cmath>
#include "SimdAlign.hpp"

namespace simdtrf {  // simdtrf namespace

auto
transform_fd(double *values, const size_t nvalues, CSimdMatrix &buffer, const size_t fd,
             const size_t nmax) -> void
{
    // NOTE: the factors are the shells' own, so they are formed once rather
    // than for every atom pair the shell pair reaches.

    const auto f_0 = 0.75 * std::sqrt(30.0);
    const auto f_1 = 0.25 * std::sqrt(30.0);
    const auto f_2 = 0.375 * std::sqrt(10.0);
    const auto f_3 = 0.75 * std::sqrt(10.0);
    const auto f_4 = 0.125 * std::sqrt(10.0);
    const auto f_5 = 0.25 * std::sqrt(10.0);
    const auto f_6 = 0.375 * std::sqrt(30.0);
    const auto f_7 = 0.125 * std::sqrt(30.0);
    const auto f_8 = 3.0 * std::sqrt(5.0);
    const auto f_9 = 0.5 * std::sqrt(15.0);
    const auto f_10 = std::sqrt(15.0);
    const auto f_11 = 1.5 * std::sqrt(5.0);
    const auto f_12 = 0.75 * std::sqrt(2.0);
    const auto f_13 = 3.0 * std::sqrt(2.0);
    const auto f_14 = 0.125 * std::sqrt(6.0);
    const auto f_15 = 0.25 * std::sqrt(6.0);
    const auto f_16 = 0.5 * std::sqrt(6.0);
    const auto f_17 = std::sqrt(6.0);
    const auto f_18 = 0.375 * std::sqrt(2.0);
    const auto f_19 = 1.5 * std::sqrt(2.0);
    const auto f_20 = 1.5 * std::sqrt(3.0);
    const auto f_21 = std::sqrt(3.0);
    const auto f_22 = 0.75 * std::sqrt(3.0);
    const auto f_23 = 0.5 * std::sqrt(3.0);
    const auto f_24 = 0.25 * std::sqrt(15.0);
    const auto f_25 = 0.75 * std::sqrt(5.0);

    // NOTE: the rows of the values are not aligned, starting at this combination's
    // offset in the values block, so they are kept out of the clause below.

    auto *g_0 = values + 0 * nvalues;
    auto *g_1 = values + 1 * nvalues;
    auto *g_2 = values + 2 * nvalues;
    auto *g_3 = values + 3 * nvalues;
    auto *g_4 = values + 4 * nvalues;
    auto *g_5 = values + 5 * nvalues;
    auto *g_6 = values + 6 * nvalues;
    auto *g_7 = values + 7 * nvalues;
    auto *g_8 = values + 8 * nvalues;
    auto *g_9 = values + 9 * nvalues;
    auto *g_10 = values + 10 * nvalues;
    auto *g_11 = values + 11 * nvalues;
    auto *g_12 = values + 12 * nvalues;
    auto *g_13 = values + 13 * nvalues;
    auto *g_14 = values + 14 * nvalues;
    auto *g_15 = values + 15 * nvalues;
    auto *g_16 = values + 16 * nvalues;
    auto *g_17 = values + 17 * nvalues;
    auto *g_18 = values + 18 * nvalues;
    auto *g_19 = values + 19 * nvalues;
    auto *g_20 = values + 20 * nvalues;
    auto *g_21 = values + 21 * nvalues;
    auto *g_22 = values + 22 * nvalues;
    auto *g_23 = values + 23 * nvalues;
    auto *g_24 = values + 24 * nvalues;
    auto *g_25 = values + 25 * nvalues;
    auto *g_26 = values + 26 * nvalues;
    auto *g_27 = values + 27 * nvalues;
    auto *g_28 = values + 28 * nvalues;
    auto *g_29 = values + 29 * nvalues;
    auto *g_30 = values + 30 * nvalues;
    auto *g_31 = values + 31 * nvalues;
    auto *g_32 = values + 32 * nvalues;
    auto *g_33 = values + 33 * nvalues;
    auto *g_34 = values + 34 * nvalues;

    const auto *fd_0 = buffer.data(fd + 0);
    const auto *fd_1 = buffer.data(fd + 1);
    const auto *fd_2 = buffer.data(fd + 2);
    const auto *fd_3 = buffer.data(fd + 3);
    const auto *fd_4 = buffer.data(fd + 4);
    const auto *fd_5 = buffer.data(fd + 5);
    const auto *fd_6 = buffer.data(fd + 6);
    const auto *fd_7 = buffer.data(fd + 7);
    const auto *fd_8 = buffer.data(fd + 8);
    const auto *fd_9 = buffer.data(fd + 9);
    const auto *fd_10 = buffer.data(fd + 10);
    const auto *fd_11 = buffer.data(fd + 11);
    const auto *fd_12 = buffer.data(fd + 12);
    const auto *fd_13 = buffer.data(fd + 13);
    const auto *fd_14 = buffer.data(fd + 14);
    const auto *fd_15 = buffer.data(fd + 15);
    const auto *fd_16 = buffer.data(fd + 16);
    const auto *fd_17 = buffer.data(fd + 17);
    const auto *fd_18 = buffer.data(fd + 18);
    const auto *fd_19 = buffer.data(fd + 19);
    const auto *fd_20 = buffer.data(fd + 20);
    const auto *fd_21 = buffer.data(fd + 21);
    const auto *fd_22 = buffer.data(fd + 22);
    const auto *fd_23 = buffer.data(fd + 23);
    const auto *fd_24 = buffer.data(fd + 24);
    const auto *fd_25 = buffer.data(fd + 25);
    const auto *fd_26 = buffer.data(fd + 26);
    const auto *fd_27 = buffer.data(fd + 27);
    const auto *fd_28 = buffer.data(fd + 28);
    const auto *fd_29 = buffer.data(fd + 29);
    const auto *fd_30 = buffer.data(fd + 30);
    const auto *fd_31 = buffer.data(fd + 31);
    const auto *fd_32 = buffer.data(fd + 32);
    const auto *fd_33 = buffer.data(fd + 33);
    const auto *fd_34 = buffer.data(fd + 34);
    const auto *fd_35 = buffer.data(fd + 35);
    const auto *fd_36 = buffer.data(fd + 36);
    const auto *fd_37 = buffer.data(fd + 37);
    const auto *fd_38 = buffer.data(fd + 38);
    const auto *fd_39 = buffer.data(fd + 39);
    const auto *fd_40 = buffer.data(fd + 40);
    const auto *fd_41 = buffer.data(fd + 41);
    const auto *fd_42 = buffer.data(fd + 42);
    const auto *fd_43 = buffer.data(fd + 43);
    const auto *fd_44 = buffer.data(fd + 44);
    const auto *fd_45 = buffer.data(fd + 45);
    const auto *fd_46 = buffer.data(fd + 46);
    const auto *fd_47 = buffer.data(fd + 47);
    const auto *fd_48 = buffer.data(fd + 48);
    const auto *fd_49 = buffer.data(fd + 49);
    const auto *fd_50 = buffer.data(fd + 50);
    const auto *fd_51 = buffer.data(fd + 51);
    const auto *fd_52 = buffer.data(fd + 52);
    const auto *fd_53 = buffer.data(fd + 53);
    const auto *fd_54 = buffer.data(fd + 54);
    const auto *fd_55 = buffer.data(fd + 55);
    const auto *fd_56 = buffer.data(fd + 56);
    const auto *fd_57 = buffer.data(fd + 57);
    const auto *fd_58 = buffer.data(fd + 58);
    const auto *fd_59 = buffer.data(fd + 59);

#pragma omp simd aligned(fd_6, fd_7, fd_8, fd_9, fd_10, fd_11, fd_36, fd_37, fd_38, fd_39, \
                         fd_40, fd_41 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_0[k] = f_0 * fd_7[k]
                 - f_1 * fd_37[k];

        g_1[k] = f_0 * fd_10[k]
                 - f_1 * fd_40[k];

        g_2[k] = -f_2 * fd_6[k]
                 - f_2 * fd_9[k]
                 + f_3 * fd_11[k]
                 + f_4 * fd_36[k]
                 + f_4 * fd_39[k]
                 - f_5 * fd_41[k];

        g_3[k] = f_0 * fd_8[k]
                 - f_1 * fd_38[k];
    }

#pragma omp simd aligned(fd_6, fd_9, fd_24, fd_25, fd_26, fd_27, fd_28, fd_29, fd_36, \
                         fd_39 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_4[k] = f_6 * fd_6[k]
                 - f_6 * fd_9[k]
                 - f_7 * fd_36[k]
                 + f_7 * fd_39[k];

        g_5[k] = f_8 * fd_25[k];

        g_6[k] = f_8 * fd_28[k];

        g_7[k] = -f_9 * fd_24[k]
                 - f_9 * fd_27[k]
                 + f_10 * fd_29[k];

        g_8[k] = f_8 * fd_26[k];

        g_9[k] = f_11 * fd_24[k]
                 - f_11 * fd_27[k];
    }

#pragma omp simd aligned(fd_7, fd_10, fd_37, fd_40, fd_49, fd_52 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_10[k] = -f_12 * fd_7[k]
                  - f_12 * fd_37[k]
                  + f_13 * fd_49[k];

        g_11[k] = -f_12 * fd_10[k]
                  - f_12 * fd_40[k]
                  + f_13 * fd_52[k];
    }

#pragma omp simd aligned(fd_6, fd_8, fd_9, fd_11, fd_36, fd_38, fd_39, fd_41, fd_48, fd_50, \
                         fd_51, fd_53 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_12[k] = f_14 * fd_6[k]
                  + f_14 * fd_9[k]
                  - f_15 * fd_11[k]
                  + f_14 * fd_36[k]
                  + f_14 * fd_39[k]
                  - f_15 * fd_41[k]
                  - f_16 * fd_48[k]
                  - f_16 * fd_51[k]
                  + f_17 * fd_53[k];

        g_13[k] = -f_12 * fd_8[k]
                  - f_12 * fd_38[k]
                  + f_13 * fd_50[k];

        g_14[k] = -f_18 * fd_6[k]
                  + f_18 * fd_9[k]
                  - f_18 * fd_36[k]
                  + f_18 * fd_39[k]
                  + f_19 * fd_48[k]
                  - f_19 * fd_51[k];
    }

#pragma omp simd aligned(fd_13, fd_16, fd_43, fd_46, fd_55, fd_58 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_15[k] = -f_20 * fd_13[k]
                  - f_20 * fd_43[k]
                  + f_21 * fd_55[k];

        g_16[k] = -f_20 * fd_16[k]
                  - f_20 * fd_46[k]
                  + f_21 * fd_58[k];
    }

#pragma omp simd aligned(fd_12, fd_14, fd_15, fd_17, fd_42, fd_44, fd_45, fd_47, fd_54, fd_56, \
                         fd_57, fd_59 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_17[k] = 0.75 * fd_12[k]
                  + 0.75 * fd_15[k]
                  - 1.5 * fd_17[k]
                  + 0.75 * fd_42[k]
                  + 0.75 * fd_45[k]
                  - 1.5 * fd_47[k]
                  - 0.5 * fd_54[k]
                  - 0.5 * fd_57[k]
                  + fd_59[k];

        g_18[k] = -f_20 * fd_14[k]
                  - f_20 * fd_44[k]
                  + f_21 * fd_56[k];

        g_19[k] = -f_22 * fd_12[k]
                  + f_22 * fd_15[k]
                  - f_22 * fd_42[k]
                  + f_22 * fd_45[k]
                  + f_23 * fd_54[k]
                  - f_23 * fd_57[k];
    }

#pragma omp simd aligned(fd_1, fd_4, fd_19, fd_22, fd_31, fd_34 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_20[k] = -f_12 * fd_1[k]
                  - f_12 * fd_19[k]
                  + f_13 * fd_31[k];

        g_21[k] = -f_12 * fd_4[k]
                  - f_12 * fd_22[k]
                  + f_13 * fd_34[k];
    }

#pragma omp simd aligned(fd_0, fd_2, fd_3, fd_5, fd_18, fd_20, fd_21, fd_23, fd_30, fd_32, \
                         fd_33, fd_35 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_22[k] = f_14 * fd_0[k]
                  + f_14 * fd_3[k]
                  - f_15 * fd_5[k]
                  + f_14 * fd_18[k]
                  + f_14 * fd_21[k]
                  - f_15 * fd_23[k]
                  - f_16 * fd_30[k]
                  - f_16 * fd_33[k]
                  + f_17 * fd_35[k];

        g_23[k] = -f_12 * fd_2[k]
                  - f_12 * fd_20[k]
                  + f_13 * fd_32[k];

        g_24[k] = -f_18 * fd_0[k]
                  + f_18 * fd_3[k]
                  - f_18 * fd_18[k]
                  + f_18 * fd_21[k]
                  + f_19 * fd_30[k]
                  - f_19 * fd_33[k];
    }

#pragma omp simd aligned(fd_12, fd_13, fd_14, fd_15, fd_16, fd_17, fd_42, fd_43, fd_44, fd_45, \
                         fd_46, fd_47 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_25[k] = f_11 * fd_13[k]
                  - f_11 * fd_43[k];

        g_26[k] = f_11 * fd_16[k]
                  - f_11 * fd_46[k];

        g_27[k] = -f_24 * fd_12[k]
                  - f_24 * fd_15[k]
                  + f_9 * fd_17[k]
                  + f_24 * fd_42[k]
                  + f_24 * fd_45[k]
                  - f_9 * fd_47[k];

        g_28[k] = f_11 * fd_14[k]
                  - f_11 * fd_44[k];
    }

#pragma omp simd aligned(fd_1, fd_4, fd_12, fd_15, fd_19, fd_22, fd_42, \
                         fd_45 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_29[k] = f_25 * fd_12[k]
                  - f_25 * fd_15[k]
                  - f_25 * fd_42[k]
                  + f_25 * fd_45[k];

        g_30[k] = f_1 * fd_1[k]
                  - f_0 * fd_19[k];

        g_31[k] = f_1 * fd_4[k]
                  - f_0 * fd_22[k];
    }

#pragma omp simd aligned(fd_0, fd_2, fd_3, fd_5, fd_18, fd_20, fd_21, \
                         fd_23 : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        g_32[k] = -f_4 * fd_0[k]
                  - f_4 * fd_3[k]
                  + f_5 * fd_5[k]
                  + f_2 * fd_18[k]
                  + f_2 * fd_21[k]
                  - f_3 * fd_23[k];

        g_33[k] = f_1 * fd_2[k]
                  - f_0 * fd_20[k];

        g_34[k] = f_7 * fd_0[k]
                  - f_7 * fd_3[k]
                  - f_6 * fd_18[k]
                  + f_6 * fd_21[k];
    }
}

}  // namespace simdtrf
