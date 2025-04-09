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

#include "TwoCenterElectronRepulsionPrimRecDD.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_dd(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_dd,
                                const size_t idx_eri_0_sd,
                                const size_t idx_eri_1_sd,
                                const size_t idx_eri_1_pp,
                                const size_t idx_eri_1_pd,
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

    // Set up components of auxiliary buffer : SD

    auto g_0_xx_0 = pbuffer.data(idx_eri_0_sd);

    auto g_0_xy_0 = pbuffer.data(idx_eri_0_sd + 1);

    auto g_0_xz_0 = pbuffer.data(idx_eri_0_sd + 2);

    auto g_0_yy_0 = pbuffer.data(idx_eri_0_sd + 3);

    auto g_0_yz_0 = pbuffer.data(idx_eri_0_sd + 4);

    auto g_0_zz_0 = pbuffer.data(idx_eri_0_sd + 5);

    // Set up components of auxiliary buffer : SD

    auto g_0_xx_1 = pbuffer.data(idx_eri_1_sd);

    auto g_0_xy_1 = pbuffer.data(idx_eri_1_sd + 1);

    auto g_0_xz_1 = pbuffer.data(idx_eri_1_sd + 2);

    auto g_0_yy_1 = pbuffer.data(idx_eri_1_sd + 3);

    auto g_0_yz_1 = pbuffer.data(idx_eri_1_sd + 4);

    auto g_0_zz_1 = pbuffer.data(idx_eri_1_sd + 5);

    // Set up components of auxiliary buffer : PP

    auto g_x_x_1 = pbuffer.data(idx_eri_1_pp);

    auto g_x_y_1 = pbuffer.data(idx_eri_1_pp + 1);

    auto g_x_z_1 = pbuffer.data(idx_eri_1_pp + 2);

    auto g_y_x_1 = pbuffer.data(idx_eri_1_pp + 3);

    auto g_y_y_1 = pbuffer.data(idx_eri_1_pp + 4);

    auto g_y_z_1 = pbuffer.data(idx_eri_1_pp + 5);

    auto g_z_x_1 = pbuffer.data(idx_eri_1_pp + 6);

    auto g_z_y_1 = pbuffer.data(idx_eri_1_pp + 7);

    auto g_z_z_1 = pbuffer.data(idx_eri_1_pp + 8);

    // Set up components of auxiliary buffer : PD

    auto g_x_xx_1 = pbuffer.data(idx_eri_1_pd);

    auto g_x_xy_1 = pbuffer.data(idx_eri_1_pd + 1);

    auto g_x_xz_1 = pbuffer.data(idx_eri_1_pd + 2);

    auto g_x_yy_1 = pbuffer.data(idx_eri_1_pd + 3);

    auto g_x_yz_1 = pbuffer.data(idx_eri_1_pd + 4);

    auto g_x_zz_1 = pbuffer.data(idx_eri_1_pd + 5);

    auto g_y_xx_1 = pbuffer.data(idx_eri_1_pd + 6);

    auto g_y_xy_1 = pbuffer.data(idx_eri_1_pd + 7);

    auto g_y_xz_1 = pbuffer.data(idx_eri_1_pd + 8);

    auto g_y_yy_1 = pbuffer.data(idx_eri_1_pd + 9);

    auto g_y_yz_1 = pbuffer.data(idx_eri_1_pd + 10);

    auto g_y_zz_1 = pbuffer.data(idx_eri_1_pd + 11);

    auto g_z_xx_1 = pbuffer.data(idx_eri_1_pd + 12);

    auto g_z_xy_1 = pbuffer.data(idx_eri_1_pd + 13);

    auto g_z_xz_1 = pbuffer.data(idx_eri_1_pd + 14);

    auto g_z_yy_1 = pbuffer.data(idx_eri_1_pd + 15);

    auto g_z_yz_1 = pbuffer.data(idx_eri_1_pd + 16);

    auto g_z_zz_1 = pbuffer.data(idx_eri_1_pd + 17);

    // Set up 0-6 components of targeted buffer : DD

    auto g_xx_xx_0 = pbuffer.data(idx_eri_0_dd);

    auto g_xx_xy_0 = pbuffer.data(idx_eri_0_dd + 1);

    auto g_xx_xz_0 = pbuffer.data(idx_eri_0_dd + 2);

    auto g_xx_yy_0 = pbuffer.data(idx_eri_0_dd + 3);

    auto g_xx_yz_0 = pbuffer.data(idx_eri_0_dd + 4);

    auto g_xx_zz_0 = pbuffer.data(idx_eri_0_dd + 5);

    #pragma omp simd aligned(g_0_xx_0, g_0_xx_1, g_0_xy_0, g_0_xy_1, g_0_xz_0, g_0_xz_1, g_0_yy_0, g_0_yy_1, g_0_yz_0, g_0_yz_1, g_0_zz_0, g_0_zz_1, g_x_x_1, g_x_xx_1, g_x_xy_1, g_x_xz_1, g_x_y_1, g_x_yy_1, g_x_yz_1, g_x_z_1, g_x_zz_1, g_xx_xx_0, g_xx_xy_0, g_xx_xz_0, g_xx_yy_0, g_xx_yz_0, g_xx_zz_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xx_xx_0[i] = g_0_xx_0[i] * fbe_0 - g_0_xx_1[i] * fz_be_0 + 2.0 * g_x_x_1[i] * fe_0 + g_x_xx_1[i] * pa_x[i];

        g_xx_xy_0[i] = g_0_xy_0[i] * fbe_0 - g_0_xy_1[i] * fz_be_0 + g_x_y_1[i] * fe_0 + g_x_xy_1[i] * pa_x[i];

        g_xx_xz_0[i] = g_0_xz_0[i] * fbe_0 - g_0_xz_1[i] * fz_be_0 + g_x_z_1[i] * fe_0 + g_x_xz_1[i] * pa_x[i];

        g_xx_yy_0[i] = g_0_yy_0[i] * fbe_0 - g_0_yy_1[i] * fz_be_0 + g_x_yy_1[i] * pa_x[i];

        g_xx_yz_0[i] = g_0_yz_0[i] * fbe_0 - g_0_yz_1[i] * fz_be_0 + g_x_yz_1[i] * pa_x[i];

        g_xx_zz_0[i] = g_0_zz_0[i] * fbe_0 - g_0_zz_1[i] * fz_be_0 + g_x_zz_1[i] * pa_x[i];
    }

    // Set up 6-12 components of targeted buffer : DD

    auto g_xy_xx_0 = pbuffer.data(idx_eri_0_dd + 6);

    auto g_xy_xy_0 = pbuffer.data(idx_eri_0_dd + 7);

    auto g_xy_xz_0 = pbuffer.data(idx_eri_0_dd + 8);

    auto g_xy_yy_0 = pbuffer.data(idx_eri_0_dd + 9);

    auto g_xy_yz_0 = pbuffer.data(idx_eri_0_dd + 10);

    auto g_xy_zz_0 = pbuffer.data(idx_eri_0_dd + 11);

    #pragma omp simd aligned(g_x_xx_1, g_x_xz_1, g_xy_xx_0, g_xy_xy_0, g_xy_xz_0, g_xy_yy_0, g_xy_yz_0, g_xy_zz_0, g_y_xy_1, g_y_y_1, g_y_yy_1, g_y_yz_1, g_y_zz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xy_xx_0[i] = g_x_xx_1[i] * pa_y[i];

        g_xy_xy_0[i] = g_y_y_1[i] * fe_0 + g_y_xy_1[i] * pa_x[i];

        g_xy_xz_0[i] = g_x_xz_1[i] * pa_y[i];

        g_xy_yy_0[i] = g_y_yy_1[i] * pa_x[i];

        g_xy_yz_0[i] = g_y_yz_1[i] * pa_x[i];

        g_xy_zz_0[i] = g_y_zz_1[i] * pa_x[i];
    }

    // Set up 12-18 components of targeted buffer : DD

    auto g_xz_xx_0 = pbuffer.data(idx_eri_0_dd + 12);

    auto g_xz_xy_0 = pbuffer.data(idx_eri_0_dd + 13);

    auto g_xz_xz_0 = pbuffer.data(idx_eri_0_dd + 14);

    auto g_xz_yy_0 = pbuffer.data(idx_eri_0_dd + 15);

    auto g_xz_yz_0 = pbuffer.data(idx_eri_0_dd + 16);

    auto g_xz_zz_0 = pbuffer.data(idx_eri_0_dd + 17);

    #pragma omp simd aligned(g_x_xx_1, g_x_xy_1, g_xz_xx_0, g_xz_xy_0, g_xz_xz_0, g_xz_yy_0, g_xz_yz_0, g_xz_zz_0, g_z_xz_1, g_z_yy_1, g_z_yz_1, g_z_z_1, g_z_zz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xz_xx_0[i] = g_x_xx_1[i] * pa_z[i];

        g_xz_xy_0[i] = g_x_xy_1[i] * pa_z[i];

        g_xz_xz_0[i] = g_z_z_1[i] * fe_0 + g_z_xz_1[i] * pa_x[i];

        g_xz_yy_0[i] = g_z_yy_1[i] * pa_x[i];

        g_xz_yz_0[i] = g_z_yz_1[i] * pa_x[i];

        g_xz_zz_0[i] = g_z_zz_1[i] * pa_x[i];
    }

    // Set up 18-24 components of targeted buffer : DD

    auto g_yy_xx_0 = pbuffer.data(idx_eri_0_dd + 18);

    auto g_yy_xy_0 = pbuffer.data(idx_eri_0_dd + 19);

    auto g_yy_xz_0 = pbuffer.data(idx_eri_0_dd + 20);

    auto g_yy_yy_0 = pbuffer.data(idx_eri_0_dd + 21);

    auto g_yy_yz_0 = pbuffer.data(idx_eri_0_dd + 22);

    auto g_yy_zz_0 = pbuffer.data(idx_eri_0_dd + 23);

    #pragma omp simd aligned(g_0_xx_0, g_0_xx_1, g_0_xy_0, g_0_xy_1, g_0_xz_0, g_0_xz_1, g_0_yy_0, g_0_yy_1, g_0_yz_0, g_0_yz_1, g_0_zz_0, g_0_zz_1, g_y_x_1, g_y_xx_1, g_y_xy_1, g_y_xz_1, g_y_y_1, g_y_yy_1, g_y_yz_1, g_y_z_1, g_y_zz_1, g_yy_xx_0, g_yy_xy_0, g_yy_xz_0, g_yy_yy_0, g_yy_yz_0, g_yy_zz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yy_xx_0[i] = g_0_xx_0[i] * fbe_0 - g_0_xx_1[i] * fz_be_0 + g_y_xx_1[i] * pa_y[i];

        g_yy_xy_0[i] = g_0_xy_0[i] * fbe_0 - g_0_xy_1[i] * fz_be_0 + g_y_x_1[i] * fe_0 + g_y_xy_1[i] * pa_y[i];

        g_yy_xz_0[i] = g_0_xz_0[i] * fbe_0 - g_0_xz_1[i] * fz_be_0 + g_y_xz_1[i] * pa_y[i];

        g_yy_yy_0[i] = g_0_yy_0[i] * fbe_0 - g_0_yy_1[i] * fz_be_0 + 2.0 * g_y_y_1[i] * fe_0 + g_y_yy_1[i] * pa_y[i];

        g_yy_yz_0[i] = g_0_yz_0[i] * fbe_0 - g_0_yz_1[i] * fz_be_0 + g_y_z_1[i] * fe_0 + g_y_yz_1[i] * pa_y[i];

        g_yy_zz_0[i] = g_0_zz_0[i] * fbe_0 - g_0_zz_1[i] * fz_be_0 + g_y_zz_1[i] * pa_y[i];
    }

    // Set up 24-30 components of targeted buffer : DD

    auto g_yz_xx_0 = pbuffer.data(idx_eri_0_dd + 24);

    auto g_yz_xy_0 = pbuffer.data(idx_eri_0_dd + 25);

    auto g_yz_xz_0 = pbuffer.data(idx_eri_0_dd + 26);

    auto g_yz_yy_0 = pbuffer.data(idx_eri_0_dd + 27);

    auto g_yz_yz_0 = pbuffer.data(idx_eri_0_dd + 28);

    auto g_yz_zz_0 = pbuffer.data(idx_eri_0_dd + 29);

    #pragma omp simd aligned(g_y_xy_1, g_y_yy_1, g_yz_xx_0, g_yz_xy_0, g_yz_xz_0, g_yz_yy_0, g_yz_yz_0, g_yz_zz_0, g_z_xx_1, g_z_xz_1, g_z_yz_1, g_z_z_1, g_z_zz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yz_xx_0[i] = g_z_xx_1[i] * pa_y[i];

        g_yz_xy_0[i] = g_y_xy_1[i] * pa_z[i];

        g_yz_xz_0[i] = g_z_xz_1[i] * pa_y[i];

        g_yz_yy_0[i] = g_y_yy_1[i] * pa_z[i];

        g_yz_yz_0[i] = g_z_z_1[i] * fe_0 + g_z_yz_1[i] * pa_y[i];

        g_yz_zz_0[i] = g_z_zz_1[i] * pa_y[i];
    }

    // Set up 30-36 components of targeted buffer : DD

    auto g_zz_xx_0 = pbuffer.data(idx_eri_0_dd + 30);

    auto g_zz_xy_0 = pbuffer.data(idx_eri_0_dd + 31);

    auto g_zz_xz_0 = pbuffer.data(idx_eri_0_dd + 32);

    auto g_zz_yy_0 = pbuffer.data(idx_eri_0_dd + 33);

    auto g_zz_yz_0 = pbuffer.data(idx_eri_0_dd + 34);

    auto g_zz_zz_0 = pbuffer.data(idx_eri_0_dd + 35);

    #pragma omp simd aligned(g_0_xx_0, g_0_xx_1, g_0_xy_0, g_0_xy_1, g_0_xz_0, g_0_xz_1, g_0_yy_0, g_0_yy_1, g_0_yz_0, g_0_yz_1, g_0_zz_0, g_0_zz_1, g_z_x_1, g_z_xx_1, g_z_xy_1, g_z_xz_1, g_z_y_1, g_z_yy_1, g_z_yz_1, g_z_z_1, g_z_zz_1, g_zz_xx_0, g_zz_xy_0, g_zz_xz_0, g_zz_yy_0, g_zz_yz_0, g_zz_zz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zz_xx_0[i] = g_0_xx_0[i] * fbe_0 - g_0_xx_1[i] * fz_be_0 + g_z_xx_1[i] * pa_z[i];

        g_zz_xy_0[i] = g_0_xy_0[i] * fbe_0 - g_0_xy_1[i] * fz_be_0 + g_z_xy_1[i] * pa_z[i];

        g_zz_xz_0[i] = g_0_xz_0[i] * fbe_0 - g_0_xz_1[i] * fz_be_0 + g_z_x_1[i] * fe_0 + g_z_xz_1[i] * pa_z[i];

        g_zz_yy_0[i] = g_0_yy_0[i] * fbe_0 - g_0_yy_1[i] * fz_be_0 + g_z_yy_1[i] * pa_z[i];

        g_zz_yz_0[i] = g_0_yz_0[i] * fbe_0 - g_0_yz_1[i] * fz_be_0 + g_z_y_1[i] * fe_0 + g_z_yz_1[i] * pa_z[i];

        g_zz_zz_0[i] = g_0_zz_0[i] * fbe_0 - g_0_zz_1[i] * fz_be_0 + 2.0 * g_z_z_1[i] * fe_0 + g_z_zz_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

