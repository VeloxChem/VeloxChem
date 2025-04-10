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

#include "TwoCenterElectronRepulsionPrimRecDP.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_dp(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_dp,
                                const size_t idx_eri_0_sp,
                                const size_t idx_eri_1_sp,
                                const size_t idx_eri_1_ps,
                                const size_t idx_eri_1_pp,
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

    // Set up components of auxiliary buffer : SP

    auto g_0_x_0 = pbuffer.data(idx_eri_0_sp);

    auto g_0_y_0 = pbuffer.data(idx_eri_0_sp + 1);

    auto g_0_z_0 = pbuffer.data(idx_eri_0_sp + 2);

    // Set up components of auxiliary buffer : SP

    auto g_0_x_1 = pbuffer.data(idx_eri_1_sp);

    auto g_0_y_1 = pbuffer.data(idx_eri_1_sp + 1);

    auto g_0_z_1 = pbuffer.data(idx_eri_1_sp + 2);

    // Set up components of auxiliary buffer : PS

    auto g_x_0_1 = pbuffer.data(idx_eri_1_ps);

    auto g_y_0_1 = pbuffer.data(idx_eri_1_ps + 1);

    auto g_z_0_1 = pbuffer.data(idx_eri_1_ps + 2);

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

    // Set up 0-3 components of targeted buffer : DP

    auto g_xx_x_0 = pbuffer.data(idx_eri_0_dp);

    auto g_xx_y_0 = pbuffer.data(idx_eri_0_dp + 1);

    auto g_xx_z_0 = pbuffer.data(idx_eri_0_dp + 2);

    #pragma omp simd aligned(g_0_x_0, g_0_x_1, g_0_y_0, g_0_y_1, g_0_z_0, g_0_z_1, g_x_0_1, g_x_x_1, g_x_y_1, g_x_z_1, g_xx_x_0, g_xx_y_0, g_xx_z_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xx_x_0[i] = g_0_x_0[i] * fbe_0 - g_0_x_1[i] * fz_be_0 + g_x_0_1[i] * fe_0 + g_x_x_1[i] * pa_x[i];

        g_xx_y_0[i] = g_0_y_0[i] * fbe_0 - g_0_y_1[i] * fz_be_0 + g_x_y_1[i] * pa_x[i];

        g_xx_z_0[i] = g_0_z_0[i] * fbe_0 - g_0_z_1[i] * fz_be_0 + g_x_z_1[i] * pa_x[i];
    }

    // Set up 3-6 components of targeted buffer : DP

    auto g_xy_x_0 = pbuffer.data(idx_eri_0_dp + 3);

    auto g_xy_y_0 = pbuffer.data(idx_eri_0_dp + 4);

    auto g_xy_z_0 = pbuffer.data(idx_eri_0_dp + 5);

    #pragma omp simd aligned(g_x_x_1, g_xy_x_0, g_xy_y_0, g_xy_z_0, g_y_y_1, g_y_z_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_xy_x_0[i] = g_x_x_1[i] * pa_y[i];

        g_xy_y_0[i] = g_y_y_1[i] * pa_x[i];

        g_xy_z_0[i] = g_y_z_1[i] * pa_x[i];
    }

    // Set up 6-9 components of targeted buffer : DP

    auto g_xz_x_0 = pbuffer.data(idx_eri_0_dp + 6);

    auto g_xz_y_0 = pbuffer.data(idx_eri_0_dp + 7);

    auto g_xz_z_0 = pbuffer.data(idx_eri_0_dp + 8);

    #pragma omp simd aligned(g_x_x_1, g_xz_x_0, g_xz_y_0, g_xz_z_0, g_z_y_1, g_z_z_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_xz_x_0[i] = g_x_x_1[i] * pa_z[i];

        g_xz_y_0[i] = g_z_y_1[i] * pa_x[i];

        g_xz_z_0[i] = g_z_z_1[i] * pa_x[i];
    }

    // Set up 9-12 components of targeted buffer : DP

    auto g_yy_x_0 = pbuffer.data(idx_eri_0_dp + 9);

    auto g_yy_y_0 = pbuffer.data(idx_eri_0_dp + 10);

    auto g_yy_z_0 = pbuffer.data(idx_eri_0_dp + 11);

    #pragma omp simd aligned(g_0_x_0, g_0_x_1, g_0_y_0, g_0_y_1, g_0_z_0, g_0_z_1, g_y_0_1, g_y_x_1, g_y_y_1, g_y_z_1, g_yy_x_0, g_yy_y_0, g_yy_z_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yy_x_0[i] = g_0_x_0[i] * fbe_0 - g_0_x_1[i] * fz_be_0 + g_y_x_1[i] * pa_y[i];

        g_yy_y_0[i] = g_0_y_0[i] * fbe_0 - g_0_y_1[i] * fz_be_0 + g_y_0_1[i] * fe_0 + g_y_y_1[i] * pa_y[i];

        g_yy_z_0[i] = g_0_z_0[i] * fbe_0 - g_0_z_1[i] * fz_be_0 + g_y_z_1[i] * pa_y[i];
    }

    // Set up 12-15 components of targeted buffer : DP

    auto g_yz_x_0 = pbuffer.data(idx_eri_0_dp + 12);

    auto g_yz_y_0 = pbuffer.data(idx_eri_0_dp + 13);

    auto g_yz_z_0 = pbuffer.data(idx_eri_0_dp + 14);

    #pragma omp simd aligned(g_y_y_1, g_yz_x_0, g_yz_y_0, g_yz_z_0, g_z_x_1, g_z_z_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_yz_x_0[i] = g_z_x_1[i] * pa_y[i];

        g_yz_y_0[i] = g_y_y_1[i] * pa_z[i];

        g_yz_z_0[i] = g_z_z_1[i] * pa_y[i];
    }

    // Set up 15-18 components of targeted buffer : DP

    auto g_zz_x_0 = pbuffer.data(idx_eri_0_dp + 15);

    auto g_zz_y_0 = pbuffer.data(idx_eri_0_dp + 16);

    auto g_zz_z_0 = pbuffer.data(idx_eri_0_dp + 17);

    #pragma omp simd aligned(g_0_x_0, g_0_x_1, g_0_y_0, g_0_y_1, g_0_z_0, g_0_z_1, g_z_0_1, g_z_x_1, g_z_y_1, g_z_z_1, g_zz_x_0, g_zz_y_0, g_zz_z_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zz_x_0[i] = g_0_x_0[i] * fbe_0 - g_0_x_1[i] * fz_be_0 + g_z_x_1[i] * pa_z[i];

        g_zz_y_0[i] = g_0_y_0[i] * fbe_0 - g_0_y_1[i] * fz_be_0 + g_z_y_1[i] * pa_z[i];

        g_zz_z_0[i] = g_0_z_0[i] * fbe_0 - g_0_z_1[i] * fz_be_0 + g_z_0_1[i] * fe_0 + g_z_z_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

