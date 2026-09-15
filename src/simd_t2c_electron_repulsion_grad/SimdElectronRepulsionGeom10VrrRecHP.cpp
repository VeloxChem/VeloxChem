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


#include "SimdElectronRepulsionGeom10VrrRecHP.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_prim_geom_10_hp_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                             const size_t gp, const size_t ip,
                                             const size_t ncols, const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_1 = buffer.data(gp + 1);
    const auto *gp_2 = buffer.data(gp + 2);
    const auto *gp_3 = buffer.data(gp + 3);
    const auto *gp_4 = buffer.data(gp + 4);
    const auto *gp_5 = buffer.data(gp + 5);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_7 = buffer.data(gp + 7);
    const auto *gp_8 = buffer.data(gp + 8);
    const auto *gp_9 = buffer.data(gp + 9);
    const auto *gp_10 = buffer.data(gp + 10);
    const auto *gp_11 = buffer.data(gp + 11);
    const auto *gp_12 = buffer.data(gp + 12);
    const auto *gp_13 = buffer.data(gp + 13);
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_15 = buffer.data(gp + 15);
    const auto *gp_16 = buffer.data(gp + 16);
    const auto *gp_17 = buffer.data(gp + 17);
    const auto *gp_18 = buffer.data(gp + 18);
    const auto *gp_19 = buffer.data(gp + 19);
    const auto *gp_20 = buffer.data(gp + 20);
    const auto *gp_21 = buffer.data(gp + 21);
    const auto *gp_22 = buffer.data(gp + 22);
    const auto *gp_23 = buffer.data(gp + 23);
    const auto *gp_24 = buffer.data(gp + 24);
    const auto *gp_25 = buffer.data(gp + 25);
    const auto *gp_26 = buffer.data(gp + 26);
    const auto *gp_27 = buffer.data(gp + 27);
    const auto *gp_28 = buffer.data(gp + 28);
    const auto *gp_29 = buffer.data(gp + 29);
    const auto *gp_30 = buffer.data(gp + 30);
    const auto *gp_31 = buffer.data(gp + 31);
    const auto *gp_32 = buffer.data(gp + 32);
    const auto *gp_33 = buffer.data(gp + 33);
    const auto *gp_34 = buffer.data(gp + 34);
    const auto *gp_35 = buffer.data(gp + 35);
    const auto *gp_36 = buffer.data(gp + 36);
    const auto *gp_37 = buffer.data(gp + 37);
    const auto *gp_38 = buffer.data(gp + 38);
    const auto *gp_39 = buffer.data(gp + 39);
    const auto *gp_40 = buffer.data(gp + 40);
    const auto *gp_41 = buffer.data(gp + 41);
    const auto *gp_42 = buffer.data(gp + 42);
    const auto *gp_43 = buffer.data(gp + 43);
    const auto *gp_44 = buffer.data(gp + 44);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_1 = buffer.data(ip + 1);
    const auto *ip_2 = buffer.data(ip + 2);
    const auto *ip_3 = buffer.data(ip + 3);
    const auto *ip_4 = buffer.data(ip + 4);
    const auto *ip_5 = buffer.data(ip + 5);
    const auto *ip_6 = buffer.data(ip + 6);
    const auto *ip_7 = buffer.data(ip + 7);
    const auto *ip_8 = buffer.data(ip + 8);
    const auto *ip_9 = buffer.data(ip + 9);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_11 = buffer.data(ip + 11);
    const auto *ip_12 = buffer.data(ip + 12);
    const auto *ip_13 = buffer.data(ip + 13);
    const auto *ip_14 = buffer.data(ip + 14);
    const auto *ip_15 = buffer.data(ip + 15);
    const auto *ip_16 = buffer.data(ip + 16);
    const auto *ip_17 = buffer.data(ip + 17);
    const auto *ip_18 = buffer.data(ip + 18);
    const auto *ip_19 = buffer.data(ip + 19);
    const auto *ip_20 = buffer.data(ip + 20);
    const auto *ip_21 = buffer.data(ip + 21);
    const auto *ip_22 = buffer.data(ip + 22);
    const auto *ip_23 = buffer.data(ip + 23);
    const auto *ip_24 = buffer.data(ip + 24);
    const auto *ip_25 = buffer.data(ip + 25);
    const auto *ip_26 = buffer.data(ip + 26);
    const auto *ip_27 = buffer.data(ip + 27);
    const auto *ip_28 = buffer.data(ip + 28);
    const auto *ip_29 = buffer.data(ip + 29);
    const auto *ip_30 = buffer.data(ip + 30);
    const auto *ip_31 = buffer.data(ip + 31);
    const auto *ip_32 = buffer.data(ip + 32);
    const auto *ip_33 = buffer.data(ip + 33);
    const auto *ip_34 = buffer.data(ip + 34);
    const auto *ip_35 = buffer.data(ip + 35);
    const auto *ip_36 = buffer.data(ip + 36);
    const auto *ip_37 = buffer.data(ip + 37);
    const auto *ip_38 = buffer.data(ip + 38);
    const auto *ip_39 = buffer.data(ip + 39);
    const auto *ip_40 = buffer.data(ip + 40);
    const auto *ip_41 = buffer.data(ip + 41);
    const auto *ip_42 = buffer.data(ip + 42);
    const auto *ip_43 = buffer.data(ip + 43);
    const auto *ip_44 = buffer.data(ip + 44);
    const auto *ip_45 = buffer.data(ip + 45);
    const auto *ip_46 = buffer.data(ip + 46);
    const auto *ip_47 = buffer.data(ip + 47);
    const auto *ip_48 = buffer.data(ip + 48);
    const auto *ip_49 = buffer.data(ip + 49);
    const auto *ip_50 = buffer.data(ip + 50);
    const auto *ip_51 = buffer.data(ip + 51);
    const auto *ip_52 = buffer.data(ip + 52);
    const auto *ip_53 = buffer.data(ip + 53);
    const auto *ip_54 = buffer.data(ip + 54);
    const auto *ip_55 = buffer.data(ip + 55);
    const auto *ip_56 = buffer.data(ip + 56);
    const auto *ip_57 = buffer.data(ip + 57);
    const auto *ip_58 = buffer.data(ip + 58);
    const auto *ip_59 = buffer.data(ip + 59);
    const auto *ip_60 = buffer.data(ip + 60);
    const auto *ip_61 = buffer.data(ip + 61);
    const auto *ip_62 = buffer.data(ip + 62);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, gp_0, gp_1, gp_2, gp_3, gp_4, ip_0, ip_1, \
                         ip_2, ip_3, ip_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = -5.0 * gp_0[k]
                 + f_0 * ip_0[k];

        t_1[k] = -5.0 * gp_1[k]
                 + f_0 * ip_1[k];

        t_2[k] = -5.0 * gp_2[k]
                 + f_0 * ip_2[k];

        t_3[k] = -4.0 * gp_3[k]
                 + f_0 * ip_3[k];

        t_4[k] = -4.0 * gp_4[k]
                 + f_0 * ip_4[k];
    }

#pragma omp simd aligned(t_5, t_6, t_7, t_8, t_9, gp_5, gp_6, gp_7, gp_8, gp_9, ip_5, ip_6, \
                         ip_7, ip_8, ip_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = -4.0 * gp_5[k]
                 + f_0 * ip_5[k];

        t_6[k] = -4.0 * gp_6[k]
                 + f_0 * ip_6[k];

        t_7[k] = -4.0 * gp_7[k]
                 + f_0 * ip_7[k];

        t_8[k] = -4.0 * gp_8[k]
                 + f_0 * ip_8[k];

        t_9[k] = -3.0 * gp_9[k]
                 + f_0 * ip_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, gp_10, gp_11, gp_12, gp_13, gp_14, \
                         ip_10, ip_11, ip_12, ip_13, ip_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -3.0 * gp_10[k]
                  + f_0 * ip_10[k];

        t_11[k] = -3.0 * gp_11[k]
                  + f_0 * ip_11[k];

        t_12[k] = -3.0 * gp_12[k]
                  + f_0 * ip_12[k];

        t_13[k] = -3.0 * gp_13[k]
                  + f_0 * ip_13[k];

        t_14[k] = -3.0 * gp_14[k]
                  + f_0 * ip_14[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, t_19, gp_15, gp_16, gp_17, gp_18, gp_19, \
                         ip_15, ip_16, ip_17, ip_18, ip_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = -3.0 * gp_15[k]
                  + f_0 * ip_15[k];

        t_16[k] = -3.0 * gp_16[k]
                  + f_0 * ip_16[k];

        t_17[k] = -3.0 * gp_17[k]
                  + f_0 * ip_17[k];

        t_18[k] = -2.0 * gp_18[k]
                  + f_0 * ip_18[k];

        t_19[k] = -2.0 * gp_19[k]
                  + f_0 * ip_19[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, t_24, gp_20, gp_21, gp_22, gp_23, gp_24, \
                         ip_20, ip_21, ip_22, ip_23, ip_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -2.0 * gp_20[k]
                  + f_0 * ip_20[k];

        t_21[k] = -2.0 * gp_21[k]
                  + f_0 * ip_21[k];

        t_22[k] = -2.0 * gp_22[k]
                  + f_0 * ip_22[k];

        t_23[k] = -2.0 * gp_23[k]
                  + f_0 * ip_23[k];

        t_24[k] = -2.0 * gp_24[k]
                  + f_0 * ip_24[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, gp_25, gp_26, gp_27, gp_28, gp_29, \
                         ip_25, ip_26, ip_27, ip_28, ip_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = -2.0 * gp_25[k]
                  + f_0 * ip_25[k];

        t_26[k] = -2.0 * gp_26[k]
                  + f_0 * ip_26[k];

        t_27[k] = -2.0 * gp_27[k]
                  + f_0 * ip_27[k];

        t_28[k] = -2.0 * gp_28[k]
                  + f_0 * ip_28[k];

        t_29[k] = -2.0 * gp_29[k]
                  + f_0 * ip_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, gp_30, gp_31, gp_32, gp_33, gp_34, \
                         ip_30, ip_31, ip_32, ip_33, ip_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -gp_30[k]
                  + f_0 * ip_30[k];

        t_31[k] = -gp_31[k]
                  + f_0 * ip_31[k];

        t_32[k] = -gp_32[k]
                  + f_0 * ip_32[k];

        t_33[k] = -gp_33[k]
                  + f_0 * ip_33[k];

        t_34[k] = -gp_34[k]
                  + f_0 * ip_34[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, gp_35, gp_36, gp_37, gp_38, gp_39, \
                         ip_35, ip_36, ip_37, ip_38, ip_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -gp_35[k]
                  + f_0 * ip_35[k];

        t_36[k] = -gp_36[k]
                  + f_0 * ip_36[k];

        t_37[k] = -gp_37[k]
                  + f_0 * ip_37[k];

        t_38[k] = -gp_38[k]
                  + f_0 * ip_38[k];

        t_39[k] = -gp_39[k]
                  + f_0 * ip_39[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, gp_40, gp_41, gp_42, gp_43, gp_44, \
                         ip_40, ip_41, ip_42, ip_43, ip_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -gp_40[k]
                  + f_0 * ip_40[k];

        t_41[k] = -gp_41[k]
                  + f_0 * ip_41[k];

        t_42[k] = -gp_42[k]
                  + f_0 * ip_42[k];

        t_43[k] = -gp_43[k]
                  + f_0 * ip_43[k];

        t_44[k] = -gp_44[k]
                  + f_0 * ip_44[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, t_50, t_51, t_52, ip_45, ip_46, ip_47, \
                         ip_48, ip_49, ip_50, ip_51, ip_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_0 * ip_45[k];

        t_46[k] = f_0 * ip_46[k];

        t_47[k] = f_0 * ip_47[k];

        t_48[k] = f_0 * ip_48[k];

        t_49[k] = f_0 * ip_49[k];

        t_50[k] = f_0 * ip_50[k];

        t_51[k] = f_0 * ip_51[k];

        t_52[k] = f_0 * ip_52[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, t_57, t_58, t_59, t_60, ip_53, ip_54, ip_55, \
                         ip_56, ip_57, ip_58, ip_59, ip_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_0 * ip_53[k];

        t_54[k] = f_0 * ip_54[k];

        t_55[k] = f_0 * ip_55[k];

        t_56[k] = f_0 * ip_56[k];

        t_57[k] = f_0 * ip_57[k];

        t_58[k] = f_0 * ip_58[k];

        t_59[k] = f_0 * ip_59[k];

        t_60[k] = f_0 * ip_60[k];
    }

#pragma omp simd aligned(t_61, t_62, ip_61, ip_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_0 * ip_61[k];

        t_62[k] = f_0 * ip_62[k];
    }
}

auto
compute_prim_geom_10_hp_electron_repulsion_1(CSimdMatrix &buffer, const size_t target,
                                             const size_t gp, const size_t ip,
                                             const size_t ncols, const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_1 = buffer.data(gp + 1);
    const auto *gp_2 = buffer.data(gp + 2);
    const auto *gp_3 = buffer.data(gp + 3);
    const auto *gp_4 = buffer.data(gp + 4);
    const auto *gp_5 = buffer.data(gp + 5);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_7 = buffer.data(gp + 7);
    const auto *gp_8 = buffer.data(gp + 8);
    const auto *gp_9 = buffer.data(gp + 9);
    const auto *gp_10 = buffer.data(gp + 10);
    const auto *gp_11 = buffer.data(gp + 11);
    const auto *gp_12 = buffer.data(gp + 12);
    const auto *gp_13 = buffer.data(gp + 13);
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_15 = buffer.data(gp + 15);
    const auto *gp_16 = buffer.data(gp + 16);
    const auto *gp_17 = buffer.data(gp + 17);
    const auto *gp_18 = buffer.data(gp + 18);
    const auto *gp_19 = buffer.data(gp + 19);
    const auto *gp_20 = buffer.data(gp + 20);
    const auto *gp_21 = buffer.data(gp + 21);
    const auto *gp_22 = buffer.data(gp + 22);
    const auto *gp_23 = buffer.data(gp + 23);
    const auto *gp_24 = buffer.data(gp + 24);
    const auto *gp_25 = buffer.data(gp + 25);
    const auto *gp_26 = buffer.data(gp + 26);
    const auto *gp_27 = buffer.data(gp + 27);
    const auto *gp_28 = buffer.data(gp + 28);
    const auto *gp_29 = buffer.data(gp + 29);
    const auto *gp_30 = buffer.data(gp + 30);
    const auto *gp_31 = buffer.data(gp + 31);
    const auto *gp_32 = buffer.data(gp + 32);
    const auto *gp_33 = buffer.data(gp + 33);
    const auto *gp_34 = buffer.data(gp + 34);
    const auto *gp_35 = buffer.data(gp + 35);
    const auto *gp_36 = buffer.data(gp + 36);
    const auto *gp_37 = buffer.data(gp + 37);
    const auto *gp_38 = buffer.data(gp + 38);
    const auto *gp_39 = buffer.data(gp + 39);
    const auto *gp_40 = buffer.data(gp + 40);
    const auto *gp_41 = buffer.data(gp + 41);
    const auto *gp_42 = buffer.data(gp + 42);
    const auto *gp_43 = buffer.data(gp + 43);
    const auto *gp_44 = buffer.data(gp + 44);

    const auto *ip_3 = buffer.data(ip + 3);
    const auto *ip_4 = buffer.data(ip + 4);
    const auto *ip_5 = buffer.data(ip + 5);
    const auto *ip_9 = buffer.data(ip + 9);
    const auto *ip_10 = buffer.data(ip + 10);
    const auto *ip_11 = buffer.data(ip + 11);
    const auto *ip_12 = buffer.data(ip + 12);
    const auto *ip_13 = buffer.data(ip + 13);
    const auto *ip_14 = buffer.data(ip + 14);
    const auto *ip_18 = buffer.data(ip + 18);
    const auto *ip_19 = buffer.data(ip + 19);
    const auto *ip_20 = buffer.data(ip + 20);
    const auto *ip_21 = buffer.data(ip + 21);
    const auto *ip_22 = buffer.data(ip + 22);
    const auto *ip_23 = buffer.data(ip + 23);
    const auto *ip_24 = buffer.data(ip + 24);
    const auto *ip_25 = buffer.data(ip + 25);
    const auto *ip_26 = buffer.data(ip + 26);
    const auto *ip_30 = buffer.data(ip + 30);
    const auto *ip_31 = buffer.data(ip + 31);
    const auto *ip_32 = buffer.data(ip + 32);
    const auto *ip_33 = buffer.data(ip + 33);
    const auto *ip_34 = buffer.data(ip + 34);
    const auto *ip_35 = buffer.data(ip + 35);
    const auto *ip_36 = buffer.data(ip + 36);
    const auto *ip_37 = buffer.data(ip + 37);
    const auto *ip_38 = buffer.data(ip + 38);
    const auto *ip_39 = buffer.data(ip + 39);
    const auto *ip_40 = buffer.data(ip + 40);
    const auto *ip_41 = buffer.data(ip + 41);
    const auto *ip_45 = buffer.data(ip + 45);
    const auto *ip_46 = buffer.data(ip + 46);
    const auto *ip_47 = buffer.data(ip + 47);
    const auto *ip_48 = buffer.data(ip + 48);
    const auto *ip_49 = buffer.data(ip + 49);
    const auto *ip_50 = buffer.data(ip + 50);
    const auto *ip_51 = buffer.data(ip + 51);
    const auto *ip_52 = buffer.data(ip + 52);
    const auto *ip_53 = buffer.data(ip + 53);
    const auto *ip_54 = buffer.data(ip + 54);
    const auto *ip_55 = buffer.data(ip + 55);
    const auto *ip_56 = buffer.data(ip + 56);
    const auto *ip_57 = buffer.data(ip + 57);
    const auto *ip_58 = buffer.data(ip + 58);
    const auto *ip_59 = buffer.data(ip + 59);
    const auto *ip_63 = buffer.data(ip + 63);
    const auto *ip_64 = buffer.data(ip + 64);
    const auto *ip_65 = buffer.data(ip + 65);
    const auto *ip_66 = buffer.data(ip + 66);
    const auto *ip_67 = buffer.data(ip + 67);
    const auto *ip_68 = buffer.data(ip + 68);
    const auto *ip_69 = buffer.data(ip + 69);
    const auto *ip_70 = buffer.data(ip + 70);
    const auto *ip_71 = buffer.data(ip + 71);
    const auto *ip_72 = buffer.data(ip + 72);
    const auto *ip_73 = buffer.data(ip + 73);
    const auto *ip_74 = buffer.data(ip + 74);
    const auto *ip_75 = buffer.data(ip + 75);
    const auto *ip_76 = buffer.data(ip + 76);
    const auto *ip_77 = buffer.data(ip + 77);
    const auto *ip_78 = buffer.data(ip + 78);
    const auto *ip_79 = buffer.data(ip + 79);
    const auto *ip_80 = buffer.data(ip + 80);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, gp_0, gp_1, gp_2, ip_3, ip_4, ip_5, \
                         ip_9, ip_10, ip_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ip_3[k];

        t_1[k] = f_0 * ip_4[k];

        t_2[k] = f_0 * ip_5[k];

        t_3[k] = -gp_0[k]
                 + f_0 * ip_9[k];

        t_4[k] = -gp_1[k]
                 + f_0 * ip_10[k];

        t_5[k] = -gp_2[k]
                 + f_0 * ip_11[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, gp_3, gp_4, gp_5, ip_12, ip_13, \
                         ip_14, ip_18, ip_19, ip_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * ip_12[k];

        t_7[k] = f_0 * ip_13[k];

        t_8[k] = f_0 * ip_14[k];

        t_9[k] = -2.0 * gp_3[k]
                 + f_0 * ip_18[k];

        t_10[k] = -2.0 * gp_4[k]
                  + f_0 * ip_19[k];

        t_11[k] = -2.0 * gp_5[k]
                  + f_0 * ip_20[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, gp_6, gp_7, gp_8, ip_21, ip_22, \
                         ip_23, ip_24, ip_25, ip_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = -gp_6[k]
                  + f_0 * ip_21[k];

        t_13[k] = -gp_7[k]
                  + f_0 * ip_22[k];

        t_14[k] = -gp_8[k]
                  + f_0 * ip_23[k];

        t_15[k] = f_0 * ip_24[k];

        t_16[k] = f_0 * ip_25[k];

        t_17[k] = f_0 * ip_26[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, gp_9, gp_10, gp_11, gp_12, gp_13, \
                         ip_30, ip_31, ip_32, ip_33, ip_34 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = -3.0 * gp_9[k]
                  + f_0 * ip_30[k];

        t_19[k] = -3.0 * gp_10[k]
                  + f_0 * ip_31[k];

        t_20[k] = -3.0 * gp_11[k]
                  + f_0 * ip_32[k];

        t_21[k] = -2.0 * gp_12[k]
                  + f_0 * ip_33[k];

        t_22[k] = -2.0 * gp_13[k]
                  + f_0 * ip_34[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, t_27, t_28, gp_14, gp_15, gp_16, gp_17, \
                         ip_35, ip_36, ip_37, ip_38, ip_39, ip_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = -2.0 * gp_14[k]
                  + f_0 * ip_35[k];

        t_24[k] = -gp_15[k]
                  + f_0 * ip_36[k];

        t_25[k] = -gp_16[k]
                  + f_0 * ip_37[k];

        t_26[k] = -gp_17[k]
                  + f_0 * ip_38[k];

        t_27[k] = f_0 * ip_39[k];

        t_28[k] = f_0 * ip_40[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, gp_18, gp_19, gp_20, gp_21, ip_41, \
                         ip_45, ip_46, ip_47, ip_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = f_0 * ip_41[k];

        t_30[k] = -4.0 * gp_18[k]
                  + f_0 * ip_45[k];

        t_31[k] = -4.0 * gp_19[k]
                  + f_0 * ip_46[k];

        t_32[k] = -4.0 * gp_20[k]
                  + f_0 * ip_47[k];

        t_33[k] = -3.0 * gp_21[k]
                  + f_0 * ip_48[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, t_37, t_38, gp_22, gp_23, gp_24, gp_25, gp_26, \
                         ip_49, ip_50, ip_51, ip_52, ip_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = -3.0 * gp_22[k]
                  + f_0 * ip_49[k];

        t_35[k] = -3.0 * gp_23[k]
                  + f_0 * ip_50[k];

        t_36[k] = -2.0 * gp_24[k]
                  + f_0 * ip_51[k];

        t_37[k] = -2.0 * gp_25[k]
                  + f_0 * ip_52[k];

        t_38[k] = -2.0 * gp_26[k]
                  + f_0 * ip_53[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, t_43, t_44, gp_27, gp_28, gp_29, ip_54, \
                         ip_55, ip_56, ip_57, ip_58, ip_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = -gp_27[k]
                  + f_0 * ip_54[k];

        t_40[k] = -gp_28[k]
                  + f_0 * ip_55[k];

        t_41[k] = -gp_29[k]
                  + f_0 * ip_56[k];

        t_42[k] = f_0 * ip_57[k];

        t_43[k] = f_0 * ip_58[k];

        t_44[k] = f_0 * ip_59[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, gp_30, gp_31, gp_32, gp_33, gp_34, \
                         ip_63, ip_64, ip_65, ip_66, ip_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = -5.0 * gp_30[k]
                  + f_0 * ip_63[k];

        t_46[k] = -5.0 * gp_31[k]
                  + f_0 * ip_64[k];

        t_47[k] = -5.0 * gp_32[k]
                  + f_0 * ip_65[k];

        t_48[k] = -4.0 * gp_33[k]
                  + f_0 * ip_66[k];

        t_49[k] = -4.0 * gp_34[k]
                  + f_0 * ip_67[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, t_54, gp_35, gp_36, gp_37, gp_38, gp_39, \
                         ip_68, ip_69, ip_70, ip_71, ip_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = -4.0 * gp_35[k]
                  + f_0 * ip_68[k];

        t_51[k] = -3.0 * gp_36[k]
                  + f_0 * ip_69[k];

        t_52[k] = -3.0 * gp_37[k]
                  + f_0 * ip_70[k];

        t_53[k] = -3.0 * gp_38[k]
                  + f_0 * ip_71[k];

        t_54[k] = -2.0 * gp_39[k]
                  + f_0 * ip_72[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, t_59, gp_40, gp_41, gp_42, gp_43, gp_44, \
                         ip_73, ip_74, ip_75, ip_76, ip_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -2.0 * gp_40[k]
                  + f_0 * ip_73[k];

        t_56[k] = -2.0 * gp_41[k]
                  + f_0 * ip_74[k];

        t_57[k] = -gp_42[k]
                  + f_0 * ip_75[k];

        t_58[k] = -gp_43[k]
                  + f_0 * ip_76[k];

        t_59[k] = -gp_44[k]
                  + f_0 * ip_77[k];
    }

#pragma omp simd aligned(t_60, t_61, t_62, ip_78, ip_79, ip_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_60[k] = f_0 * ip_78[k];

        t_61[k] = f_0 * ip_79[k];

        t_62[k] = f_0 * ip_80[k];
    }
}

auto
compute_prim_geom_10_hp_electron_repulsion_2(CSimdMatrix &buffer, const size_t target,
                                             const size_t gp, const size_t ip,
                                             const size_t ncols, const double alpha) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 * alpha;

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

    const auto *gp_0 = buffer.data(gp + 0);
    const auto *gp_1 = buffer.data(gp + 1);
    const auto *gp_2 = buffer.data(gp + 2);
    const auto *gp_3 = buffer.data(gp + 3);
    const auto *gp_4 = buffer.data(gp + 4);
    const auto *gp_5 = buffer.data(gp + 5);
    const auto *gp_6 = buffer.data(gp + 6);
    const auto *gp_7 = buffer.data(gp + 7);
    const auto *gp_8 = buffer.data(gp + 8);
    const auto *gp_9 = buffer.data(gp + 9);
    const auto *gp_10 = buffer.data(gp + 10);
    const auto *gp_11 = buffer.data(gp + 11);
    const auto *gp_12 = buffer.data(gp + 12);
    const auto *gp_13 = buffer.data(gp + 13);
    const auto *gp_14 = buffer.data(gp + 14);
    const auto *gp_15 = buffer.data(gp + 15);
    const auto *gp_16 = buffer.data(gp + 16);
    const auto *gp_17 = buffer.data(gp + 17);
    const auto *gp_18 = buffer.data(gp + 18);
    const auto *gp_19 = buffer.data(gp + 19);
    const auto *gp_20 = buffer.data(gp + 20);
    const auto *gp_21 = buffer.data(gp + 21);
    const auto *gp_22 = buffer.data(gp + 22);
    const auto *gp_23 = buffer.data(gp + 23);
    const auto *gp_24 = buffer.data(gp + 24);
    const auto *gp_25 = buffer.data(gp + 25);
    const auto *gp_26 = buffer.data(gp + 26);
    const auto *gp_27 = buffer.data(gp + 27);
    const auto *gp_28 = buffer.data(gp + 28);
    const auto *gp_29 = buffer.data(gp + 29);
    const auto *gp_30 = buffer.data(gp + 30);
    const auto *gp_31 = buffer.data(gp + 31);
    const auto *gp_32 = buffer.data(gp + 32);
    const auto *gp_33 = buffer.data(gp + 33);
    const auto *gp_34 = buffer.data(gp + 34);
    const auto *gp_35 = buffer.data(gp + 35);
    const auto *gp_36 = buffer.data(gp + 36);
    const auto *gp_37 = buffer.data(gp + 37);
    const auto *gp_38 = buffer.data(gp + 38);
    const auto *gp_39 = buffer.data(gp + 39);
    const auto *gp_40 = buffer.data(gp + 40);
    const auto *gp_41 = buffer.data(gp + 41);
    const auto *gp_42 = buffer.data(gp + 42);
    const auto *gp_43 = buffer.data(gp + 43);
    const auto *gp_44 = buffer.data(gp + 44);

    const auto *ip_6 = buffer.data(ip + 6);
    const auto *ip_7 = buffer.data(ip + 7);
    const auto *ip_8 = buffer.data(ip + 8);
    const auto *ip_12 = buffer.data(ip + 12);
    const auto *ip_13 = buffer.data(ip + 13);
    const auto *ip_14 = buffer.data(ip + 14);
    const auto *ip_15 = buffer.data(ip + 15);
    const auto *ip_16 = buffer.data(ip + 16);
    const auto *ip_17 = buffer.data(ip + 17);
    const auto *ip_21 = buffer.data(ip + 21);
    const auto *ip_22 = buffer.data(ip + 22);
    const auto *ip_23 = buffer.data(ip + 23);
    const auto *ip_24 = buffer.data(ip + 24);
    const auto *ip_25 = buffer.data(ip + 25);
    const auto *ip_26 = buffer.data(ip + 26);
    const auto *ip_27 = buffer.data(ip + 27);
    const auto *ip_28 = buffer.data(ip + 28);
    const auto *ip_29 = buffer.data(ip + 29);
    const auto *ip_33 = buffer.data(ip + 33);
    const auto *ip_34 = buffer.data(ip + 34);
    const auto *ip_35 = buffer.data(ip + 35);
    const auto *ip_36 = buffer.data(ip + 36);
    const auto *ip_37 = buffer.data(ip + 37);
    const auto *ip_38 = buffer.data(ip + 38);
    const auto *ip_39 = buffer.data(ip + 39);
    const auto *ip_40 = buffer.data(ip + 40);
    const auto *ip_41 = buffer.data(ip + 41);
    const auto *ip_42 = buffer.data(ip + 42);
    const auto *ip_43 = buffer.data(ip + 43);
    const auto *ip_44 = buffer.data(ip + 44);
    const auto *ip_48 = buffer.data(ip + 48);
    const auto *ip_49 = buffer.data(ip + 49);
    const auto *ip_50 = buffer.data(ip + 50);
    const auto *ip_51 = buffer.data(ip + 51);
    const auto *ip_52 = buffer.data(ip + 52);
    const auto *ip_53 = buffer.data(ip + 53);
    const auto *ip_54 = buffer.data(ip + 54);
    const auto *ip_55 = buffer.data(ip + 55);
    const auto *ip_56 = buffer.data(ip + 56);
    const auto *ip_57 = buffer.data(ip + 57);
    const auto *ip_58 = buffer.data(ip + 58);
    const auto *ip_59 = buffer.data(ip + 59);
    const auto *ip_60 = buffer.data(ip + 60);
    const auto *ip_61 = buffer.data(ip + 61);
    const auto *ip_62 = buffer.data(ip + 62);
    const auto *ip_66 = buffer.data(ip + 66);
    const auto *ip_67 = buffer.data(ip + 67);
    const auto *ip_68 = buffer.data(ip + 68);
    const auto *ip_69 = buffer.data(ip + 69);
    const auto *ip_70 = buffer.data(ip + 70);
    const auto *ip_71 = buffer.data(ip + 71);
    const auto *ip_72 = buffer.data(ip + 72);
    const auto *ip_73 = buffer.data(ip + 73);
    const auto *ip_74 = buffer.data(ip + 74);
    const auto *ip_75 = buffer.data(ip + 75);
    const auto *ip_76 = buffer.data(ip + 76);
    const auto *ip_77 = buffer.data(ip + 77);
    const auto *ip_78 = buffer.data(ip + 78);
    const auto *ip_79 = buffer.data(ip + 79);
    const auto *ip_80 = buffer.data(ip + 80);
    const auto *ip_81 = buffer.data(ip + 81);
    const auto *ip_82 = buffer.data(ip + 82);
    const auto *ip_83 = buffer.data(ip + 83);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, t_6, gp_0, ip_6, ip_7, ip_8, ip_12, \
                         ip_13, ip_14, ip_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ip_6[k];

        t_1[k] = f_0 * ip_7[k];

        t_2[k] = f_0 * ip_8[k];

        t_3[k] = f_0 * ip_12[k];

        t_4[k] = f_0 * ip_13[k];

        t_5[k] = f_0 * ip_14[k];

        t_6[k] = -gp_0[k]
                 + f_0 * ip_15[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, t_10, t_11, t_12, gp_1, gp_2, gp_3, ip_16, ip_17, \
                         ip_21, ip_22, ip_23, ip_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = -gp_1[k]
                 + f_0 * ip_16[k];

        t_8[k] = -gp_2[k]
                 + f_0 * ip_17[k];

        t_9[k] = f_0 * ip_21[k];

        t_10[k] = f_0 * ip_22[k];

        t_11[k] = f_0 * ip_23[k];

        t_12[k] = -gp_3[k]
                  + f_0 * ip_24[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, t_16, t_17, gp_4, gp_5, gp_6, gp_7, gp_8, ip_25, \
                         ip_26, ip_27, ip_28, ip_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = -gp_4[k]
                  + f_0 * ip_25[k];

        t_14[k] = -gp_5[k]
                  + f_0 * ip_26[k];

        t_15[k] = -2.0 * gp_6[k]
                  + f_0 * ip_27[k];

        t_16[k] = -2.0 * gp_7[k]
                  + f_0 * ip_28[k];

        t_17[k] = -2.0 * gp_8[k]
                  + f_0 * ip_29[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, t_21, t_22, t_23, gp_9, gp_10, gp_11, ip_33, ip_34, \
                         ip_35, ip_36, ip_37, ip_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_0 * ip_33[k];

        t_19[k] = f_0 * ip_34[k];

        t_20[k] = f_0 * ip_35[k];

        t_21[k] = -gp_9[k]
                  + f_0 * ip_36[k];

        t_22[k] = -gp_10[k]
                  + f_0 * ip_37[k];

        t_23[k] = -gp_11[k]
                  + f_0 * ip_38[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, t_28, gp_12, gp_13, gp_14, gp_15, gp_16, \
                         ip_39, ip_40, ip_41, ip_42, ip_43 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = -2.0 * gp_12[k]
                  + f_0 * ip_39[k];

        t_25[k] = -2.0 * gp_13[k]
                  + f_0 * ip_40[k];

        t_26[k] = -2.0 * gp_14[k]
                  + f_0 * ip_41[k];

        t_27[k] = -3.0 * gp_15[k]
                  + f_0 * ip_42[k];

        t_28[k] = -3.0 * gp_16[k]
                  + f_0 * ip_43[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, t_33, t_34, gp_17, gp_18, gp_19, ip_44, \
                         ip_48, ip_49, ip_50, ip_51, ip_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = -3.0 * gp_17[k]
                  + f_0 * ip_44[k];

        t_30[k] = f_0 * ip_48[k];

        t_31[k] = f_0 * ip_49[k];

        t_32[k] = f_0 * ip_50[k];

        t_33[k] = -gp_18[k]
                  + f_0 * ip_51[k];

        t_34[k] = -gp_19[k]
                  + f_0 * ip_52[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, t_39, gp_20, gp_21, gp_22, gp_23, gp_24, \
                         ip_53, ip_54, ip_55, ip_56, ip_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = -gp_20[k]
                  + f_0 * ip_53[k];

        t_36[k] = -2.0 * gp_21[k]
                  + f_0 * ip_54[k];

        t_37[k] = -2.0 * gp_22[k]
                  + f_0 * ip_55[k];

        t_38[k] = -2.0 * gp_23[k]
                  + f_0 * ip_56[k];

        t_39[k] = -3.0 * gp_24[k]
                  + f_0 * ip_57[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, t_44, gp_25, gp_26, gp_27, gp_28, gp_29, \
                         ip_58, ip_59, ip_60, ip_61, ip_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = -3.0 * gp_25[k]
                  + f_0 * ip_58[k];

        t_41[k] = -3.0 * gp_26[k]
                  + f_0 * ip_59[k];

        t_42[k] = -4.0 * gp_27[k]
                  + f_0 * ip_60[k];

        t_43[k] = -4.0 * gp_28[k]
                  + f_0 * ip_61[k];

        t_44[k] = -4.0 * gp_29[k]
                  + f_0 * ip_62[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, t_49, t_50, gp_30, gp_31, gp_32, ip_66, \
                         ip_67, ip_68, ip_69, ip_70, ip_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_0 * ip_66[k];

        t_46[k] = f_0 * ip_67[k];

        t_47[k] = f_0 * ip_68[k];

        t_48[k] = -gp_30[k]
                  + f_0 * ip_69[k];

        t_49[k] = -gp_31[k]
                  + f_0 * ip_70[k];

        t_50[k] = -gp_32[k]
                  + f_0 * ip_71[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, t_55, gp_33, gp_34, gp_35, gp_36, gp_37, \
                         ip_72, ip_73, ip_74, ip_75, ip_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = -2.0 * gp_33[k]
                  + f_0 * ip_72[k];

        t_52[k] = -2.0 * gp_34[k]
                  + f_0 * ip_73[k];

        t_53[k] = -2.0 * gp_35[k]
                  + f_0 * ip_74[k];

        t_54[k] = -3.0 * gp_36[k]
                  + f_0 * ip_75[k];

        t_55[k] = -3.0 * gp_37[k]
                  + f_0 * ip_76[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, gp_38, gp_39, gp_40, gp_41, gp_42, \
                         ip_77, ip_78, ip_79, ip_80, ip_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = -3.0 * gp_38[k]
                  + f_0 * ip_77[k];

        t_57[k] = -4.0 * gp_39[k]
                  + f_0 * ip_78[k];

        t_58[k] = -4.0 * gp_40[k]
                  + f_0 * ip_79[k];

        t_59[k] = -4.0 * gp_41[k]
                  + f_0 * ip_80[k];

        t_60[k] = -5.0 * gp_42[k]
                  + f_0 * ip_81[k];
    }

#pragma omp simd aligned(t_61, t_62, gp_43, gp_44, ip_82, ip_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = -5.0 * gp_43[k]
                  + f_0 * ip_82[k];

        t_62[k] = -5.0 * gp_44[k]
                  + f_0 * ip_83[k];
    }
}

}  // namespace simdt2ceri
