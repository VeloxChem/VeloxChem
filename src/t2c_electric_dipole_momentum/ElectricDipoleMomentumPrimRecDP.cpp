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

#include "ElectricDipoleMomentumPrimRecDP.hpp"

namespace diprec {  // diprec namespace

auto
comp_prim_electric_dipole_momentum_dp(CSimdArray<double>&       pbuffer,
                                      const size_t              idx_dip_dp,
                                      const size_t              idx_dip_sp,
                                      const size_t              idx_dip_ps,
                                      const size_t              idx_ovl_pp,
                                      const size_t              idx_dip_pp,
                                      const CSimdArray<double>& factors,
                                      const size_t              idx_rpa,
                                      const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : SP

    auto tr_x_0_x = pbuffer.data(idx_dip_sp);

    auto tr_x_0_y = pbuffer.data(idx_dip_sp + 1);

    auto tr_x_0_z = pbuffer.data(idx_dip_sp + 2);

    auto tr_y_0_x = pbuffer.data(idx_dip_sp + 3);

    auto tr_y_0_y = pbuffer.data(idx_dip_sp + 4);

    auto tr_y_0_z = pbuffer.data(idx_dip_sp + 5);

    auto tr_z_0_x = pbuffer.data(idx_dip_sp + 6);

    auto tr_z_0_y = pbuffer.data(idx_dip_sp + 7);

    auto tr_z_0_z = pbuffer.data(idx_dip_sp + 8);

    // Set up components of auxiliary buffer : PS

    auto tr_x_x_0 = pbuffer.data(idx_dip_ps);

    auto tr_x_y_0 = pbuffer.data(idx_dip_ps + 1);

    auto tr_x_z_0 = pbuffer.data(idx_dip_ps + 2);

    auto tr_y_x_0 = pbuffer.data(idx_dip_ps + 3);

    auto tr_y_y_0 = pbuffer.data(idx_dip_ps + 4);

    auto tr_y_z_0 = pbuffer.data(idx_dip_ps + 5);

    auto tr_z_x_0 = pbuffer.data(idx_dip_ps + 6);

    auto tr_z_y_0 = pbuffer.data(idx_dip_ps + 7);

    auto tr_z_z_0 = pbuffer.data(idx_dip_ps + 8);

    // Set up components of auxiliary buffer : PP

    auto ts_x_x = pbuffer.data(idx_ovl_pp);

    auto ts_x_y = pbuffer.data(idx_ovl_pp + 1);

    auto ts_x_z = pbuffer.data(idx_ovl_pp + 2);

    auto ts_y_x = pbuffer.data(idx_ovl_pp + 3);

    auto ts_y_y = pbuffer.data(idx_ovl_pp + 4);

    auto ts_y_z = pbuffer.data(idx_ovl_pp + 5);

    auto ts_z_x = pbuffer.data(idx_ovl_pp + 6);

    auto ts_z_y = pbuffer.data(idx_ovl_pp + 7);

    auto ts_z_z = pbuffer.data(idx_ovl_pp + 8);

    // Set up components of auxiliary buffer : PP

    auto tr_x_x_x = pbuffer.data(idx_dip_pp);

    auto tr_x_x_y = pbuffer.data(idx_dip_pp + 1);

    auto tr_x_x_z = pbuffer.data(idx_dip_pp + 2);

    auto tr_x_y_x = pbuffer.data(idx_dip_pp + 3);

    auto tr_x_y_y = pbuffer.data(idx_dip_pp + 4);

    auto tr_x_y_z = pbuffer.data(idx_dip_pp + 5);

    auto tr_x_z_x = pbuffer.data(idx_dip_pp + 6);

    auto tr_x_z_y = pbuffer.data(idx_dip_pp + 7);

    auto tr_x_z_z = pbuffer.data(idx_dip_pp + 8);

    auto tr_y_x_x = pbuffer.data(idx_dip_pp + 9);

    auto tr_y_x_y = pbuffer.data(idx_dip_pp + 10);

    auto tr_y_x_z = pbuffer.data(idx_dip_pp + 11);

    auto tr_y_y_x = pbuffer.data(idx_dip_pp + 12);

    auto tr_y_y_y = pbuffer.data(idx_dip_pp + 13);

    auto tr_y_y_z = pbuffer.data(idx_dip_pp + 14);

    auto tr_y_z_x = pbuffer.data(idx_dip_pp + 15);

    auto tr_y_z_y = pbuffer.data(idx_dip_pp + 16);

    auto tr_y_z_z = pbuffer.data(idx_dip_pp + 17);

    auto tr_z_x_x = pbuffer.data(idx_dip_pp + 18);

    auto tr_z_x_y = pbuffer.data(idx_dip_pp + 19);

    auto tr_z_x_z = pbuffer.data(idx_dip_pp + 20);

    auto tr_z_y_x = pbuffer.data(idx_dip_pp + 21);

    auto tr_z_y_y = pbuffer.data(idx_dip_pp + 22);

    auto tr_z_y_z = pbuffer.data(idx_dip_pp + 23);

    auto tr_z_z_x = pbuffer.data(idx_dip_pp + 24);

    auto tr_z_z_y = pbuffer.data(idx_dip_pp + 25);

    auto tr_z_z_z = pbuffer.data(idx_dip_pp + 26);

    // Set up 0-3 components of targeted buffer : DP

    auto tr_x_xx_x = pbuffer.data(idx_dip_dp);

    auto tr_x_xx_y = pbuffer.data(idx_dip_dp + 1);

    auto tr_x_xx_z = pbuffer.data(idx_dip_dp + 2);

#pragma omp simd aligned(pa_x,          \
                             tr_x_0_x,  \
                             tr_x_0_y,  \
                             tr_x_0_z,  \
                             tr_x_x_0,  \
                             tr_x_x_x,  \
                             tr_x_x_y,  \
                             tr_x_x_z,  \
                             tr_x_xx_x, \
                             tr_x_xx_y, \
                             tr_x_xx_z, \
                             ts_x_x,    \
                             ts_x_y,    \
                             ts_x_z,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xx_x[i] = tr_x_0_x[i] * fe_0 + tr_x_x_0[i] * fe_0 + ts_x_x[i] * fe_0 + tr_x_x_x[i] * pa_x[i];

        tr_x_xx_y[i] = tr_x_0_y[i] * fe_0 + ts_x_y[i] * fe_0 + tr_x_x_y[i] * pa_x[i];

        tr_x_xx_z[i] = tr_x_0_z[i] * fe_0 + ts_x_z[i] * fe_0 + tr_x_x_z[i] * pa_x[i];
    }

    // Set up 3-6 components of targeted buffer : DP

    auto tr_x_xy_x = pbuffer.data(idx_dip_dp + 3);

    auto tr_x_xy_y = pbuffer.data(idx_dip_dp + 4);

    auto tr_x_xy_z = pbuffer.data(idx_dip_dp + 5);

#pragma omp simd aligned(pa_x, pa_y, tr_x_x_x, tr_x_x_z, tr_x_xy_x, tr_x_xy_y, tr_x_xy_z, tr_x_y_y, ts_y_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xy_x[i] = tr_x_x_x[i] * pa_y[i];

        tr_x_xy_y[i] = ts_y_y[i] * fe_0 + tr_x_y_y[i] * pa_x[i];

        tr_x_xy_z[i] = tr_x_x_z[i] * pa_y[i];
    }

    // Set up 6-9 components of targeted buffer : DP

    auto tr_x_xz_x = pbuffer.data(idx_dip_dp + 6);

    auto tr_x_xz_y = pbuffer.data(idx_dip_dp + 7);

    auto tr_x_xz_z = pbuffer.data(idx_dip_dp + 8);

#pragma omp simd aligned(pa_x, pa_z, tr_x_x_x, tr_x_x_y, tr_x_xz_x, tr_x_xz_y, tr_x_xz_z, tr_x_z_z, ts_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xz_x[i] = tr_x_x_x[i] * pa_z[i];

        tr_x_xz_y[i] = tr_x_x_y[i] * pa_z[i];

        tr_x_xz_z[i] = ts_z_z[i] * fe_0 + tr_x_z_z[i] * pa_x[i];
    }

    // Set up 9-12 components of targeted buffer : DP

    auto tr_x_yy_x = pbuffer.data(idx_dip_dp + 9);

    auto tr_x_yy_y = pbuffer.data(idx_dip_dp + 10);

    auto tr_x_yy_z = pbuffer.data(idx_dip_dp + 11);

#pragma omp simd aligned(pa_y, tr_x_0_x, tr_x_0_y, tr_x_0_z, tr_x_y_0, tr_x_y_x, tr_x_y_y, tr_x_y_z, tr_x_yy_x, tr_x_yy_y, tr_x_yy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yy_x[i] = tr_x_0_x[i] * fe_0 + tr_x_y_x[i] * pa_y[i];

        tr_x_yy_y[i] = tr_x_0_y[i] * fe_0 + tr_x_y_0[i] * fe_0 + tr_x_y_y[i] * pa_y[i];

        tr_x_yy_z[i] = tr_x_0_z[i] * fe_0 + tr_x_y_z[i] * pa_y[i];
    }

    // Set up 12-15 components of targeted buffer : DP

    auto tr_x_yz_x = pbuffer.data(idx_dip_dp + 12);

    auto tr_x_yz_y = pbuffer.data(idx_dip_dp + 13);

    auto tr_x_yz_z = pbuffer.data(idx_dip_dp + 14);

#pragma omp simd aligned(pa_y, pa_z, tr_x_y_y, tr_x_yz_x, tr_x_yz_y, tr_x_yz_z, tr_x_z_x, tr_x_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tr_x_yz_x[i] = tr_x_z_x[i] * pa_y[i];

        tr_x_yz_y[i] = tr_x_y_y[i] * pa_z[i];

        tr_x_yz_z[i] = tr_x_z_z[i] * pa_y[i];
    }

    // Set up 15-18 components of targeted buffer : DP

    auto tr_x_zz_x = pbuffer.data(idx_dip_dp + 15);

    auto tr_x_zz_y = pbuffer.data(idx_dip_dp + 16);

    auto tr_x_zz_z = pbuffer.data(idx_dip_dp + 17);

#pragma omp simd aligned(pa_z, tr_x_0_x, tr_x_0_y, tr_x_0_z, tr_x_z_0, tr_x_z_x, tr_x_z_y, tr_x_z_z, tr_x_zz_x, tr_x_zz_y, tr_x_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_zz_x[i] = tr_x_0_x[i] * fe_0 + tr_x_z_x[i] * pa_z[i];

        tr_x_zz_y[i] = tr_x_0_y[i] * fe_0 + tr_x_z_y[i] * pa_z[i];

        tr_x_zz_z[i] = tr_x_0_z[i] * fe_0 + tr_x_z_0[i] * fe_0 + tr_x_z_z[i] * pa_z[i];
    }

    // Set up 18-21 components of targeted buffer : DP

    auto tr_y_xx_x = pbuffer.data(idx_dip_dp + 18);

    auto tr_y_xx_y = pbuffer.data(idx_dip_dp + 19);

    auto tr_y_xx_z = pbuffer.data(idx_dip_dp + 20);

#pragma omp simd aligned(pa_x, tr_y_0_x, tr_y_0_y, tr_y_0_z, tr_y_x_0, tr_y_x_x, tr_y_x_y, tr_y_x_z, tr_y_xx_x, tr_y_xx_y, tr_y_xx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xx_x[i] = tr_y_0_x[i] * fe_0 + tr_y_x_0[i] * fe_0 + tr_y_x_x[i] * pa_x[i];

        tr_y_xx_y[i] = tr_y_0_y[i] * fe_0 + tr_y_x_y[i] * pa_x[i];

        tr_y_xx_z[i] = tr_y_0_z[i] * fe_0 + tr_y_x_z[i] * pa_x[i];
    }

    // Set up 21-24 components of targeted buffer : DP

    auto tr_y_xy_x = pbuffer.data(idx_dip_dp + 21);

    auto tr_y_xy_y = pbuffer.data(idx_dip_dp + 22);

    auto tr_y_xy_z = pbuffer.data(idx_dip_dp + 23);

#pragma omp simd aligned(pa_x, tr_y_xy_x, tr_y_xy_y, tr_y_xy_z, tr_y_y_0, tr_y_y_x, tr_y_y_y, tr_y_y_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xy_x[i] = tr_y_y_0[i] * fe_0 + tr_y_y_x[i] * pa_x[i];

        tr_y_xy_y[i] = tr_y_y_y[i] * pa_x[i];

        tr_y_xy_z[i] = tr_y_y_z[i] * pa_x[i];
    }

    // Set up 24-27 components of targeted buffer : DP

    auto tr_y_xz_x = pbuffer.data(idx_dip_dp + 24);

    auto tr_y_xz_y = pbuffer.data(idx_dip_dp + 25);

    auto tr_y_xz_z = pbuffer.data(idx_dip_dp + 26);

#pragma omp simd aligned(pa_x, pa_z, tr_y_x_x, tr_y_xz_x, tr_y_xz_y, tr_y_xz_z, tr_y_z_y, tr_y_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tr_y_xz_x[i] = tr_y_x_x[i] * pa_z[i];

        tr_y_xz_y[i] = tr_y_z_y[i] * pa_x[i];

        tr_y_xz_z[i] = tr_y_z_z[i] * pa_x[i];
    }

    // Set up 27-30 components of targeted buffer : DP

    auto tr_y_yy_x = pbuffer.data(idx_dip_dp + 27);

    auto tr_y_yy_y = pbuffer.data(idx_dip_dp + 28);

    auto tr_y_yy_z = pbuffer.data(idx_dip_dp + 29);

#pragma omp simd aligned(pa_y,          \
                             tr_y_0_x,  \
                             tr_y_0_y,  \
                             tr_y_0_z,  \
                             tr_y_y_0,  \
                             tr_y_y_x,  \
                             tr_y_y_y,  \
                             tr_y_y_z,  \
                             tr_y_yy_x, \
                             tr_y_yy_y, \
                             tr_y_yy_z, \
                             ts_y_x,    \
                             ts_y_y,    \
                             ts_y_z,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yy_x[i] = tr_y_0_x[i] * fe_0 + ts_y_x[i] * fe_0 + tr_y_y_x[i] * pa_y[i];

        tr_y_yy_y[i] = tr_y_0_y[i] * fe_0 + tr_y_y_0[i] * fe_0 + ts_y_y[i] * fe_0 + tr_y_y_y[i] * pa_y[i];

        tr_y_yy_z[i] = tr_y_0_z[i] * fe_0 + ts_y_z[i] * fe_0 + tr_y_y_z[i] * pa_y[i];
    }

    // Set up 30-33 components of targeted buffer : DP

    auto tr_y_yz_x = pbuffer.data(idx_dip_dp + 30);

    auto tr_y_yz_y = pbuffer.data(idx_dip_dp + 31);

    auto tr_y_yz_z = pbuffer.data(idx_dip_dp + 32);

#pragma omp simd aligned(pa_y, pa_z, tr_y_y_x, tr_y_y_y, tr_y_yz_x, tr_y_yz_y, tr_y_yz_z, tr_y_z_z, ts_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yz_x[i] = tr_y_y_x[i] * pa_z[i];

        tr_y_yz_y[i] = tr_y_y_y[i] * pa_z[i];

        tr_y_yz_z[i] = ts_z_z[i] * fe_0 + tr_y_z_z[i] * pa_y[i];
    }

    // Set up 33-36 components of targeted buffer : DP

    auto tr_y_zz_x = pbuffer.data(idx_dip_dp + 33);

    auto tr_y_zz_y = pbuffer.data(idx_dip_dp + 34);

    auto tr_y_zz_z = pbuffer.data(idx_dip_dp + 35);

#pragma omp simd aligned(pa_z, tr_y_0_x, tr_y_0_y, tr_y_0_z, tr_y_z_0, tr_y_z_x, tr_y_z_y, tr_y_z_z, tr_y_zz_x, tr_y_zz_y, tr_y_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_zz_x[i] = tr_y_0_x[i] * fe_0 + tr_y_z_x[i] * pa_z[i];

        tr_y_zz_y[i] = tr_y_0_y[i] * fe_0 + tr_y_z_y[i] * pa_z[i];

        tr_y_zz_z[i] = tr_y_0_z[i] * fe_0 + tr_y_z_0[i] * fe_0 + tr_y_z_z[i] * pa_z[i];
    }

    // Set up 36-39 components of targeted buffer : DP

    auto tr_z_xx_x = pbuffer.data(idx_dip_dp + 36);

    auto tr_z_xx_y = pbuffer.data(idx_dip_dp + 37);

    auto tr_z_xx_z = pbuffer.data(idx_dip_dp + 38);

#pragma omp simd aligned(pa_x, tr_z_0_x, tr_z_0_y, tr_z_0_z, tr_z_x_0, tr_z_x_x, tr_z_x_y, tr_z_x_z, tr_z_xx_x, tr_z_xx_y, tr_z_xx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xx_x[i] = tr_z_0_x[i] * fe_0 + tr_z_x_0[i] * fe_0 + tr_z_x_x[i] * pa_x[i];

        tr_z_xx_y[i] = tr_z_0_y[i] * fe_0 + tr_z_x_y[i] * pa_x[i];

        tr_z_xx_z[i] = tr_z_0_z[i] * fe_0 + tr_z_x_z[i] * pa_x[i];
    }

    // Set up 39-42 components of targeted buffer : DP

    auto tr_z_xy_x = pbuffer.data(idx_dip_dp + 39);

    auto tr_z_xy_y = pbuffer.data(idx_dip_dp + 40);

    auto tr_z_xy_z = pbuffer.data(idx_dip_dp + 41);

#pragma omp simd aligned(pa_x, pa_y, tr_z_x_x, tr_z_xy_x, tr_z_xy_y, tr_z_xy_z, tr_z_y_y, tr_z_y_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tr_z_xy_x[i] = tr_z_x_x[i] * pa_y[i];

        tr_z_xy_y[i] = tr_z_y_y[i] * pa_x[i];

        tr_z_xy_z[i] = tr_z_y_z[i] * pa_x[i];
    }

    // Set up 42-45 components of targeted buffer : DP

    auto tr_z_xz_x = pbuffer.data(idx_dip_dp + 42);

    auto tr_z_xz_y = pbuffer.data(idx_dip_dp + 43);

    auto tr_z_xz_z = pbuffer.data(idx_dip_dp + 44);

#pragma omp simd aligned(pa_x, tr_z_xz_x, tr_z_xz_y, tr_z_xz_z, tr_z_z_0, tr_z_z_x, tr_z_z_y, tr_z_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xz_x[i] = tr_z_z_0[i] * fe_0 + tr_z_z_x[i] * pa_x[i];

        tr_z_xz_y[i] = tr_z_z_y[i] * pa_x[i];

        tr_z_xz_z[i] = tr_z_z_z[i] * pa_x[i];
    }

    // Set up 45-48 components of targeted buffer : DP

    auto tr_z_yy_x = pbuffer.data(idx_dip_dp + 45);

    auto tr_z_yy_y = pbuffer.data(idx_dip_dp + 46);

    auto tr_z_yy_z = pbuffer.data(idx_dip_dp + 47);

#pragma omp simd aligned(pa_y, tr_z_0_x, tr_z_0_y, tr_z_0_z, tr_z_y_0, tr_z_y_x, tr_z_y_y, tr_z_y_z, tr_z_yy_x, tr_z_yy_y, tr_z_yy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yy_x[i] = tr_z_0_x[i] * fe_0 + tr_z_y_x[i] * pa_y[i];

        tr_z_yy_y[i] = tr_z_0_y[i] * fe_0 + tr_z_y_0[i] * fe_0 + tr_z_y_y[i] * pa_y[i];

        tr_z_yy_z[i] = tr_z_0_z[i] * fe_0 + tr_z_y_z[i] * pa_y[i];
    }

    // Set up 48-51 components of targeted buffer : DP

    auto tr_z_yz_x = pbuffer.data(idx_dip_dp + 48);

    auto tr_z_yz_y = pbuffer.data(idx_dip_dp + 49);

    auto tr_z_yz_z = pbuffer.data(idx_dip_dp + 50);

#pragma omp simd aligned(pa_y, tr_z_yz_x, tr_z_yz_y, tr_z_yz_z, tr_z_z_0, tr_z_z_x, tr_z_z_y, tr_z_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yz_x[i] = tr_z_z_x[i] * pa_y[i];

        tr_z_yz_y[i] = tr_z_z_0[i] * fe_0 + tr_z_z_y[i] * pa_y[i];

        tr_z_yz_z[i] = tr_z_z_z[i] * pa_y[i];
    }

    // Set up 51-54 components of targeted buffer : DP

    auto tr_z_zz_x = pbuffer.data(idx_dip_dp + 51);

    auto tr_z_zz_y = pbuffer.data(idx_dip_dp + 52);

    auto tr_z_zz_z = pbuffer.data(idx_dip_dp + 53);

#pragma omp simd aligned(pa_z,          \
                             tr_z_0_x,  \
                             tr_z_0_y,  \
                             tr_z_0_z,  \
                             tr_z_z_0,  \
                             tr_z_z_x,  \
                             tr_z_z_y,  \
                             tr_z_z_z,  \
                             tr_z_zz_x, \
                             tr_z_zz_y, \
                             tr_z_zz_z, \
                             ts_z_x,    \
                             ts_z_y,    \
                             ts_z_z,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_zz_x[i] = tr_z_0_x[i] * fe_0 + ts_z_x[i] * fe_0 + tr_z_z_x[i] * pa_z[i];

        tr_z_zz_y[i] = tr_z_0_y[i] * fe_0 + ts_z_y[i] * fe_0 + tr_z_z_y[i] * pa_z[i];

        tr_z_zz_z[i] = tr_z_0_z[i] * fe_0 + tr_z_z_0[i] * fe_0 + ts_z_z[i] * fe_0 + tr_z_z_z[i] * pa_z[i];
    }
}

}  // namespace diprec
