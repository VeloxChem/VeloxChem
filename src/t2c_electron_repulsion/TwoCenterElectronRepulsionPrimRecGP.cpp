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

#include "TwoCenterElectronRepulsionPrimRecGP.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_gp(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_gp,
                                const size_t idx_eri_0_dp,
                                const size_t idx_eri_1_dp,
                                const size_t idx_eri_1_fs,
                                const size_t idx_eri_1_fp,
                                const CSimdArray<double>& factors,
                                const size_t idx_rpa,
                                const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : DP

    auto g_xx_x_0 = pbuffer.data(idx_eri_0_dp);

    auto g_xx_y_0 = pbuffer.data(idx_eri_0_dp + 1);

    auto g_xx_z_0 = pbuffer.data(idx_eri_0_dp + 2);

    auto g_yy_x_0 = pbuffer.data(idx_eri_0_dp + 9);

    auto g_yy_y_0 = pbuffer.data(idx_eri_0_dp + 10);

    auto g_yy_z_0 = pbuffer.data(idx_eri_0_dp + 11);

    auto g_zz_x_0 = pbuffer.data(idx_eri_0_dp + 15);

    auto g_zz_y_0 = pbuffer.data(idx_eri_0_dp + 16);

    auto g_zz_z_0 = pbuffer.data(idx_eri_0_dp + 17);

    // Set up components of auxiliary buffer : DP

    auto g_xx_x_1 = pbuffer.data(idx_eri_1_dp);

    auto g_xx_y_1 = pbuffer.data(idx_eri_1_dp + 1);

    auto g_xx_z_1 = pbuffer.data(idx_eri_1_dp + 2);

    auto g_yy_x_1 = pbuffer.data(idx_eri_1_dp + 9);

    auto g_yy_y_1 = pbuffer.data(idx_eri_1_dp + 10);

    auto g_yy_z_1 = pbuffer.data(idx_eri_1_dp + 11);

    auto g_zz_x_1 = pbuffer.data(idx_eri_1_dp + 15);

    auto g_zz_y_1 = pbuffer.data(idx_eri_1_dp + 16);

    auto g_zz_z_1 = pbuffer.data(idx_eri_1_dp + 17);

    // Set up components of auxiliary buffer : FS

    auto g_xxx_0_1 = pbuffer.data(idx_eri_1_fs);

    auto g_yyy_0_1 = pbuffer.data(idx_eri_1_fs + 6);

    auto g_zzz_0_1 = pbuffer.data(idx_eri_1_fs + 9);

    // Set up components of auxiliary buffer : FP

    auto g_xxx_x_1 = pbuffer.data(idx_eri_1_fp);

    auto g_xxx_y_1 = pbuffer.data(idx_eri_1_fp + 1);

    auto g_xxx_z_1 = pbuffer.data(idx_eri_1_fp + 2);

    auto g_xxy_x_1 = pbuffer.data(idx_eri_1_fp + 3);

    auto g_xxy_y_1 = pbuffer.data(idx_eri_1_fp + 4);

    auto g_xxz_x_1 = pbuffer.data(idx_eri_1_fp + 6);

    auto g_xxz_z_1 = pbuffer.data(idx_eri_1_fp + 8);

    auto g_xyy_x_1 = pbuffer.data(idx_eri_1_fp + 9);

    auto g_xyy_y_1 = pbuffer.data(idx_eri_1_fp + 10);

    auto g_xyy_z_1 = pbuffer.data(idx_eri_1_fp + 11);

    auto g_xzz_x_1 = pbuffer.data(idx_eri_1_fp + 15);

    auto g_xzz_y_1 = pbuffer.data(idx_eri_1_fp + 16);

    auto g_xzz_z_1 = pbuffer.data(idx_eri_1_fp + 17);

    auto g_yyy_x_1 = pbuffer.data(idx_eri_1_fp + 18);

    auto g_yyy_y_1 = pbuffer.data(idx_eri_1_fp + 19);

    auto g_yyy_z_1 = pbuffer.data(idx_eri_1_fp + 20);

    auto g_yyz_y_1 = pbuffer.data(idx_eri_1_fp + 22);

    auto g_yyz_z_1 = pbuffer.data(idx_eri_1_fp + 23);

    auto g_yzz_x_1 = pbuffer.data(idx_eri_1_fp + 24);

    auto g_yzz_y_1 = pbuffer.data(idx_eri_1_fp + 25);

    auto g_yzz_z_1 = pbuffer.data(idx_eri_1_fp + 26);

    auto g_zzz_x_1 = pbuffer.data(idx_eri_1_fp + 27);

    auto g_zzz_y_1 = pbuffer.data(idx_eri_1_fp + 28);

    auto g_zzz_z_1 = pbuffer.data(idx_eri_1_fp + 29);

    // Set up 0-3 components of targeted buffer : GP

    auto g_xxxx_x_0 = pbuffer.data(idx_eri_0_gp);

    auto g_xxxx_y_0 = pbuffer.data(idx_eri_0_gp + 1);

    auto g_xxxx_z_0 = pbuffer.data(idx_eri_0_gp + 2);

    #pragma omp simd aligned(g_xx_x_0, g_xx_x_1, g_xx_y_0, g_xx_y_1, g_xx_z_0, g_xx_z_1, g_xxx_0_1, g_xxx_x_1, g_xxx_y_1, g_xxx_z_1, g_xxxx_x_0, g_xxxx_y_0, g_xxxx_z_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxx_x_0[i] = 3.0 * g_xx_x_0[i] * fbe_0 - 3.0 * g_xx_x_1[i] * fz_be_0 + g_xxx_0_1[i] * fe_0 + g_xxx_x_1[i] * pa_x[i];

        g_xxxx_y_0[i] = 3.0 * g_xx_y_0[i] * fbe_0 - 3.0 * g_xx_y_1[i] * fz_be_0 + g_xxx_y_1[i] * pa_x[i];

        g_xxxx_z_0[i] = 3.0 * g_xx_z_0[i] * fbe_0 - 3.0 * g_xx_z_1[i] * fz_be_0 + g_xxx_z_1[i] * pa_x[i];
    }

    // Set up 3-6 components of targeted buffer : GP

    auto g_xxxy_x_0 = pbuffer.data(idx_eri_0_gp + 3);

    auto g_xxxy_y_0 = pbuffer.data(idx_eri_0_gp + 4);

    auto g_xxxy_z_0 = pbuffer.data(idx_eri_0_gp + 5);

    #pragma omp simd aligned(g_xxx_0_1, g_xxx_x_1, g_xxx_y_1, g_xxx_z_1, g_xxxy_x_0, g_xxxy_y_0, g_xxxy_z_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxy_x_0[i] = g_xxx_x_1[i] * pa_y[i];

        g_xxxy_y_0[i] = g_xxx_0_1[i] * fe_0 + g_xxx_y_1[i] * pa_y[i];

        g_xxxy_z_0[i] = g_xxx_z_1[i] * pa_y[i];
    }

    // Set up 6-9 components of targeted buffer : GP

    auto g_xxxz_x_0 = pbuffer.data(idx_eri_0_gp + 6);

    auto g_xxxz_y_0 = pbuffer.data(idx_eri_0_gp + 7);

    auto g_xxxz_z_0 = pbuffer.data(idx_eri_0_gp + 8);

    #pragma omp simd aligned(g_xxx_0_1, g_xxx_x_1, g_xxx_y_1, g_xxx_z_1, g_xxxz_x_0, g_xxxz_y_0, g_xxxz_z_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxz_x_0[i] = g_xxx_x_1[i] * pa_z[i];

        g_xxxz_y_0[i] = g_xxx_y_1[i] * pa_z[i];

        g_xxxz_z_0[i] = g_xxx_0_1[i] * fe_0 + g_xxx_z_1[i] * pa_z[i];
    }

    // Set up 9-12 components of targeted buffer : GP

    auto g_xxyy_x_0 = pbuffer.data(idx_eri_0_gp + 9);

    auto g_xxyy_y_0 = pbuffer.data(idx_eri_0_gp + 10);

    auto g_xxyy_z_0 = pbuffer.data(idx_eri_0_gp + 11);

    #pragma omp simd aligned(g_xx_x_0, g_xx_x_1, g_xxy_x_1, g_xxyy_x_0, g_xxyy_y_0, g_xxyy_z_0, g_xyy_y_1, g_xyy_z_1, g_yy_y_0, g_yy_y_1, g_yy_z_0, g_yy_z_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_xxyy_x_0[i] = g_xx_x_0[i] * fbe_0 - g_xx_x_1[i] * fz_be_0 + g_xxy_x_1[i] * pa_y[i];

        g_xxyy_y_0[i] = g_yy_y_0[i] * fbe_0 - g_yy_y_1[i] * fz_be_0 + g_xyy_y_1[i] * pa_x[i];

        g_xxyy_z_0[i] = g_yy_z_0[i] * fbe_0 - g_yy_z_1[i] * fz_be_0 + g_xyy_z_1[i] * pa_x[i];
    }

    // Set up 12-15 components of targeted buffer : GP

    auto g_xxyz_x_0 = pbuffer.data(idx_eri_0_gp + 12);

    auto g_xxyz_y_0 = pbuffer.data(idx_eri_0_gp + 13);

    auto g_xxyz_z_0 = pbuffer.data(idx_eri_0_gp + 14);

    #pragma omp simd aligned(g_xxy_y_1, g_xxyz_x_0, g_xxyz_y_0, g_xxyz_z_0, g_xxz_x_1, g_xxz_z_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_xxyz_x_0[i] = g_xxz_x_1[i] * pa_y[i];

        g_xxyz_y_0[i] = g_xxy_y_1[i] * pa_z[i];

        g_xxyz_z_0[i] = g_xxz_z_1[i] * pa_y[i];
    }

    // Set up 15-18 components of targeted buffer : GP

    auto g_xxzz_x_0 = pbuffer.data(idx_eri_0_gp + 15);

    auto g_xxzz_y_0 = pbuffer.data(idx_eri_0_gp + 16);

    auto g_xxzz_z_0 = pbuffer.data(idx_eri_0_gp + 17);

    #pragma omp simd aligned(g_xx_x_0, g_xx_x_1, g_xxz_x_1, g_xxzz_x_0, g_xxzz_y_0, g_xxzz_z_0, g_xzz_y_1, g_xzz_z_1, g_zz_y_0, g_zz_y_1, g_zz_z_0, g_zz_z_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_xxzz_x_0[i] = g_xx_x_0[i] * fbe_0 - g_xx_x_1[i] * fz_be_0 + g_xxz_x_1[i] * pa_z[i];

        g_xxzz_y_0[i] = g_zz_y_0[i] * fbe_0 - g_zz_y_1[i] * fz_be_0 + g_xzz_y_1[i] * pa_x[i];

        g_xxzz_z_0[i] = g_zz_z_0[i] * fbe_0 - g_zz_z_1[i] * fz_be_0 + g_xzz_z_1[i] * pa_x[i];
    }

    // Set up 18-21 components of targeted buffer : GP

    auto g_xyyy_x_0 = pbuffer.data(idx_eri_0_gp + 18);

    auto g_xyyy_y_0 = pbuffer.data(idx_eri_0_gp + 19);

    auto g_xyyy_z_0 = pbuffer.data(idx_eri_0_gp + 20);

    #pragma omp simd aligned(g_xyyy_x_0, g_xyyy_y_0, g_xyyy_z_0, g_yyy_0_1, g_yyy_x_1, g_yyy_y_1, g_yyy_z_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyy_x_0[i] = g_yyy_0_1[i] * fe_0 + g_yyy_x_1[i] * pa_x[i];

        g_xyyy_y_0[i] = g_yyy_y_1[i] * pa_x[i];

        g_xyyy_z_0[i] = g_yyy_z_1[i] * pa_x[i];
    }

    // Set up 21-24 components of targeted buffer : GP

    auto g_xyyz_x_0 = pbuffer.data(idx_eri_0_gp + 21);

    auto g_xyyz_y_0 = pbuffer.data(idx_eri_0_gp + 22);

    auto g_xyyz_z_0 = pbuffer.data(idx_eri_0_gp + 23);

    #pragma omp simd aligned(g_xyy_x_1, g_xyyz_x_0, g_xyyz_y_0, g_xyyz_z_0, g_yyz_y_1, g_yyz_z_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_xyyz_x_0[i] = g_xyy_x_1[i] * pa_z[i];

        g_xyyz_y_0[i] = g_yyz_y_1[i] * pa_x[i];

        g_xyyz_z_0[i] = g_yyz_z_1[i] * pa_x[i];
    }

    // Set up 24-27 components of targeted buffer : GP

    auto g_xyzz_x_0 = pbuffer.data(idx_eri_0_gp + 24);

    auto g_xyzz_y_0 = pbuffer.data(idx_eri_0_gp + 25);

    auto g_xyzz_z_0 = pbuffer.data(idx_eri_0_gp + 26);

    #pragma omp simd aligned(g_xyzz_x_0, g_xyzz_y_0, g_xyzz_z_0, g_xzz_x_1, g_yzz_y_1, g_yzz_z_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_xyzz_x_0[i] = g_xzz_x_1[i] * pa_y[i];

        g_xyzz_y_0[i] = g_yzz_y_1[i] * pa_x[i];

        g_xyzz_z_0[i] = g_yzz_z_1[i] * pa_x[i];
    }

    // Set up 27-30 components of targeted buffer : GP

    auto g_xzzz_x_0 = pbuffer.data(idx_eri_0_gp + 27);

    auto g_xzzz_y_0 = pbuffer.data(idx_eri_0_gp + 28);

    auto g_xzzz_z_0 = pbuffer.data(idx_eri_0_gp + 29);

    #pragma omp simd aligned(g_xzzz_x_0, g_xzzz_y_0, g_xzzz_z_0, g_zzz_0_1, g_zzz_x_1, g_zzz_y_1, g_zzz_z_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xzzz_x_0[i] = g_zzz_0_1[i] * fe_0 + g_zzz_x_1[i] * pa_x[i];

        g_xzzz_y_0[i] = g_zzz_y_1[i] * pa_x[i];

        g_xzzz_z_0[i] = g_zzz_z_1[i] * pa_x[i];
    }

    // Set up 30-33 components of targeted buffer : GP

    auto g_yyyy_x_0 = pbuffer.data(idx_eri_0_gp + 30);

    auto g_yyyy_y_0 = pbuffer.data(idx_eri_0_gp + 31);

    auto g_yyyy_z_0 = pbuffer.data(idx_eri_0_gp + 32);

    #pragma omp simd aligned(g_yy_x_0, g_yy_x_1, g_yy_y_0, g_yy_y_1, g_yy_z_0, g_yy_z_1, g_yyy_0_1, g_yyy_x_1, g_yyy_y_1, g_yyy_z_1, g_yyyy_x_0, g_yyyy_y_0, g_yyyy_z_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyy_x_0[i] = 3.0 * g_yy_x_0[i] * fbe_0 - 3.0 * g_yy_x_1[i] * fz_be_0 + g_yyy_x_1[i] * pa_y[i];

        g_yyyy_y_0[i] = 3.0 * g_yy_y_0[i] * fbe_0 - 3.0 * g_yy_y_1[i] * fz_be_0 + g_yyy_0_1[i] * fe_0 + g_yyy_y_1[i] * pa_y[i];

        g_yyyy_z_0[i] = 3.0 * g_yy_z_0[i] * fbe_0 - 3.0 * g_yy_z_1[i] * fz_be_0 + g_yyy_z_1[i] * pa_y[i];
    }

    // Set up 33-36 components of targeted buffer : GP

    auto g_yyyz_x_0 = pbuffer.data(idx_eri_0_gp + 33);

    auto g_yyyz_y_0 = pbuffer.data(idx_eri_0_gp + 34);

    auto g_yyyz_z_0 = pbuffer.data(idx_eri_0_gp + 35);

    #pragma omp simd aligned(g_yyy_0_1, g_yyy_x_1, g_yyy_y_1, g_yyy_z_1, g_yyyz_x_0, g_yyyz_y_0, g_yyyz_z_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yyyz_x_0[i] = g_yyy_x_1[i] * pa_z[i];

        g_yyyz_y_0[i] = g_yyy_y_1[i] * pa_z[i];

        g_yyyz_z_0[i] = g_yyy_0_1[i] * fe_0 + g_yyy_z_1[i] * pa_z[i];
    }

    // Set up 36-39 components of targeted buffer : GP

    auto g_yyzz_x_0 = pbuffer.data(idx_eri_0_gp + 36);

    auto g_yyzz_y_0 = pbuffer.data(idx_eri_0_gp + 37);

    auto g_yyzz_z_0 = pbuffer.data(idx_eri_0_gp + 38);

    #pragma omp simd aligned(g_yy_y_0, g_yy_y_1, g_yyz_y_1, g_yyzz_x_0, g_yyzz_y_0, g_yyzz_z_0, g_yzz_x_1, g_yzz_z_1, g_zz_x_0, g_zz_x_1, g_zz_z_0, g_zz_z_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_yyzz_x_0[i] = g_zz_x_0[i] * fbe_0 - g_zz_x_1[i] * fz_be_0 + g_yzz_x_1[i] * pa_y[i];

        g_yyzz_y_0[i] = g_yy_y_0[i] * fbe_0 - g_yy_y_1[i] * fz_be_0 + g_yyz_y_1[i] * pa_z[i];

        g_yyzz_z_0[i] = g_zz_z_0[i] * fbe_0 - g_zz_z_1[i] * fz_be_0 + g_yzz_z_1[i] * pa_y[i];
    }

    // Set up 39-42 components of targeted buffer : GP

    auto g_yzzz_x_0 = pbuffer.data(idx_eri_0_gp + 39);

    auto g_yzzz_y_0 = pbuffer.data(idx_eri_0_gp + 40);

    auto g_yzzz_z_0 = pbuffer.data(idx_eri_0_gp + 41);

    #pragma omp simd aligned(g_yzzz_x_0, g_yzzz_y_0, g_yzzz_z_0, g_zzz_0_1, g_zzz_x_1, g_zzz_y_1, g_zzz_z_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yzzz_x_0[i] = g_zzz_x_1[i] * pa_y[i];

        g_yzzz_y_0[i] = g_zzz_0_1[i] * fe_0 + g_zzz_y_1[i] * pa_y[i];

        g_yzzz_z_0[i] = g_zzz_z_1[i] * pa_y[i];
    }

    // Set up 42-45 components of targeted buffer : GP

    auto g_zzzz_x_0 = pbuffer.data(idx_eri_0_gp + 42);

    auto g_zzzz_y_0 = pbuffer.data(idx_eri_0_gp + 43);

    auto g_zzzz_z_0 = pbuffer.data(idx_eri_0_gp + 44);

    #pragma omp simd aligned(g_zz_x_0, g_zz_x_1, g_zz_y_0, g_zz_y_1, g_zz_z_0, g_zz_z_1, g_zzz_0_1, g_zzz_x_1, g_zzz_y_1, g_zzz_z_1, g_zzzz_x_0, g_zzzz_y_0, g_zzzz_z_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zzzz_x_0[i] = 3.0 * g_zz_x_0[i] * fbe_0 - 3.0 * g_zz_x_1[i] * fz_be_0 + g_zzz_x_1[i] * pa_z[i];

        g_zzzz_y_0[i] = 3.0 * g_zz_y_0[i] * fbe_0 - 3.0 * g_zz_y_1[i] * fz_be_0 + g_zzz_y_1[i] * pa_z[i];

        g_zzzz_z_0[i] = 3.0 * g_zz_z_0[i] * fbe_0 - 3.0 * g_zz_z_1[i] * fz_be_0 + g_zzz_0_1[i] * fe_0 + g_zzz_z_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

