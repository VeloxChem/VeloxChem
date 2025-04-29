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

#include "KineticEnergyPrimRecDP.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_dp(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_dp,
                            const size_t              idx_ovl_sp,
                            const size_t              idx_kin_sp,
                            const size_t              idx_kin_ps,
                            const size_t              idx_kin_pp,
                            const size_t              idx_ovl_dp,
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

    auto ts_0_x = pbuffer.data(idx_ovl_sp);

    auto ts_0_y = pbuffer.data(idx_ovl_sp + 1);

    auto ts_0_z = pbuffer.data(idx_ovl_sp + 2);

    // Set up components of auxiliary buffer : SP

    auto tk_0_x = pbuffer.data(idx_kin_sp);

    auto tk_0_y = pbuffer.data(idx_kin_sp + 1);

    auto tk_0_z = pbuffer.data(idx_kin_sp + 2);

    // Set up components of auxiliary buffer : PS

    auto tk_x_0 = pbuffer.data(idx_kin_ps);

    auto tk_y_0 = pbuffer.data(idx_kin_ps + 1);

    auto tk_z_0 = pbuffer.data(idx_kin_ps + 2);

    // Set up components of auxiliary buffer : PP

    auto tk_x_x = pbuffer.data(idx_kin_pp);

    auto tk_x_y = pbuffer.data(idx_kin_pp + 1);

    auto tk_x_z = pbuffer.data(idx_kin_pp + 2);

    auto tk_y_x = pbuffer.data(idx_kin_pp + 3);

    auto tk_y_y = pbuffer.data(idx_kin_pp + 4);

    auto tk_y_z = pbuffer.data(idx_kin_pp + 5);

    auto tk_z_x = pbuffer.data(idx_kin_pp + 6);

    auto tk_z_y = pbuffer.data(idx_kin_pp + 7);

    auto tk_z_z = pbuffer.data(idx_kin_pp + 8);

    // Set up components of auxiliary buffer : DP

    auto ts_xx_x = pbuffer.data(idx_ovl_dp);

    auto ts_xx_y = pbuffer.data(idx_ovl_dp + 1);

    auto ts_xx_z = pbuffer.data(idx_ovl_dp + 2);

    auto ts_xy_x = pbuffer.data(idx_ovl_dp + 3);

    auto ts_xy_y = pbuffer.data(idx_ovl_dp + 4);

    auto ts_xy_z = pbuffer.data(idx_ovl_dp + 5);

    auto ts_xz_x = pbuffer.data(idx_ovl_dp + 6);

    auto ts_xz_y = pbuffer.data(idx_ovl_dp + 7);

    auto ts_xz_z = pbuffer.data(idx_ovl_dp + 8);

    auto ts_yy_x = pbuffer.data(idx_ovl_dp + 9);

    auto ts_yy_y = pbuffer.data(idx_ovl_dp + 10);

    auto ts_yy_z = pbuffer.data(idx_ovl_dp + 11);

    auto ts_yz_x = pbuffer.data(idx_ovl_dp + 12);

    auto ts_yz_y = pbuffer.data(idx_ovl_dp + 13);

    auto ts_yz_z = pbuffer.data(idx_ovl_dp + 14);

    auto ts_zz_x = pbuffer.data(idx_ovl_dp + 15);

    auto ts_zz_y = pbuffer.data(idx_ovl_dp + 16);

    auto ts_zz_z = pbuffer.data(idx_ovl_dp + 17);

    // Set up 0-3 components of targeted buffer : DP

    auto tk_xx_x = pbuffer.data(idx_kin_dp);

    auto tk_xx_y = pbuffer.data(idx_kin_dp + 1);

    auto tk_xx_z = pbuffer.data(idx_kin_dp + 2);

#pragma omp simd aligned(pa_x,        \
                             tk_0_x,  \
                             tk_0_y,  \
                             tk_0_z,  \
                             tk_x_0,  \
                             tk_x_x,  \
                             tk_x_y,  \
                             tk_x_z,  \
                             tk_xx_x, \
                             tk_xx_y, \
                             tk_xx_z, \
                             ts_0_x,  \
                             ts_0_y,  \
                             ts_0_z,  \
                             ts_xx_x, \
                             ts_xx_y, \
                             ts_xx_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xx_x[i] = -2.0 * ts_0_x[i] * fbe_0 * fz_0 + tk_0_x[i] * fe_0 + tk_x_0[i] * fe_0 + tk_x_x[i] * pa_x[i] + 2.0 * ts_xx_x[i] * fz_0;

        tk_xx_y[i] = -2.0 * ts_0_y[i] * fbe_0 * fz_0 + tk_0_y[i] * fe_0 + tk_x_y[i] * pa_x[i] + 2.0 * ts_xx_y[i] * fz_0;

        tk_xx_z[i] = -2.0 * ts_0_z[i] * fbe_0 * fz_0 + tk_0_z[i] * fe_0 + tk_x_z[i] * pa_x[i] + 2.0 * ts_xx_z[i] * fz_0;
    }

    // Set up 3-6 components of targeted buffer : DP

    auto tk_xy_x = pbuffer.data(idx_kin_dp + 3);

    auto tk_xy_y = pbuffer.data(idx_kin_dp + 4);

    auto tk_xy_z = pbuffer.data(idx_kin_dp + 5);

#pragma omp simd aligned(pa_x, pa_y, tk_x_x, tk_xy_x, tk_xy_y, tk_xy_z, tk_y_y, tk_y_z, ts_xy_x, ts_xy_y, ts_xy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fz_0 = a_exp * b_exps[i] / (a_exp + b_exps[i]);

        tk_xy_x[i] = tk_x_x[i] * pa_y[i] + 2.0 * ts_xy_x[i] * fz_0;

        tk_xy_y[i] = tk_y_y[i] * pa_x[i] + 2.0 * ts_xy_y[i] * fz_0;

        tk_xy_z[i] = tk_y_z[i] * pa_x[i] + 2.0 * ts_xy_z[i] * fz_0;
    }

    // Set up 6-9 components of targeted buffer : DP

    auto tk_xz_x = pbuffer.data(idx_kin_dp + 6);

    auto tk_xz_y = pbuffer.data(idx_kin_dp + 7);

    auto tk_xz_z = pbuffer.data(idx_kin_dp + 8);

#pragma omp simd aligned(pa_x, pa_z, tk_x_x, tk_xz_x, tk_xz_y, tk_xz_z, tk_z_y, tk_z_z, ts_xz_x, ts_xz_y, ts_xz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fz_0 = a_exp * b_exps[i] / (a_exp + b_exps[i]);

        tk_xz_x[i] = tk_x_x[i] * pa_z[i] + 2.0 * ts_xz_x[i] * fz_0;

        tk_xz_y[i] = tk_z_y[i] * pa_x[i] + 2.0 * ts_xz_y[i] * fz_0;

        tk_xz_z[i] = tk_z_z[i] * pa_x[i] + 2.0 * ts_xz_z[i] * fz_0;
    }

    // Set up 9-12 components of targeted buffer : DP

    auto tk_yy_x = pbuffer.data(idx_kin_dp + 9);

    auto tk_yy_y = pbuffer.data(idx_kin_dp + 10);

    auto tk_yy_z = pbuffer.data(idx_kin_dp + 11);

#pragma omp simd aligned(pa_y,        \
                             tk_0_x,  \
                             tk_0_y,  \
                             tk_0_z,  \
                             tk_y_0,  \
                             tk_y_x,  \
                             tk_y_y,  \
                             tk_y_z,  \
                             tk_yy_x, \
                             tk_yy_y, \
                             tk_yy_z, \
                             ts_0_x,  \
                             ts_0_y,  \
                             ts_0_z,  \
                             ts_yy_x, \
                             ts_yy_y, \
                             ts_yy_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yy_x[i] = -2.0 * ts_0_x[i] * fbe_0 * fz_0 + tk_0_x[i] * fe_0 + tk_y_x[i] * pa_y[i] + 2.0 * ts_yy_x[i] * fz_0;

        tk_yy_y[i] = -2.0 * ts_0_y[i] * fbe_0 * fz_0 + tk_0_y[i] * fe_0 + tk_y_0[i] * fe_0 + tk_y_y[i] * pa_y[i] + 2.0 * ts_yy_y[i] * fz_0;

        tk_yy_z[i] = -2.0 * ts_0_z[i] * fbe_0 * fz_0 + tk_0_z[i] * fe_0 + tk_y_z[i] * pa_y[i] + 2.0 * ts_yy_z[i] * fz_0;
    }

    // Set up 12-15 components of targeted buffer : DP

    auto tk_yz_x = pbuffer.data(idx_kin_dp + 12);

    auto tk_yz_y = pbuffer.data(idx_kin_dp + 13);

    auto tk_yz_z = pbuffer.data(idx_kin_dp + 14);

#pragma omp simd aligned(pa_y, pa_z, tk_y_y, tk_yz_x, tk_yz_y, tk_yz_z, tk_z_x, tk_z_z, ts_yz_x, ts_yz_y, ts_yz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fz_0 = a_exp * b_exps[i] / (a_exp + b_exps[i]);

        tk_yz_x[i] = tk_z_x[i] * pa_y[i] + 2.0 * ts_yz_x[i] * fz_0;

        tk_yz_y[i] = tk_y_y[i] * pa_z[i] + 2.0 * ts_yz_y[i] * fz_0;

        tk_yz_z[i] = tk_z_z[i] * pa_y[i] + 2.0 * ts_yz_z[i] * fz_0;
    }

    // Set up 15-18 components of targeted buffer : DP

    auto tk_zz_x = pbuffer.data(idx_kin_dp + 15);

    auto tk_zz_y = pbuffer.data(idx_kin_dp + 16);

    auto tk_zz_z = pbuffer.data(idx_kin_dp + 17);

#pragma omp simd aligned(pa_z,        \
                             tk_0_x,  \
                             tk_0_y,  \
                             tk_0_z,  \
                             tk_z_0,  \
                             tk_z_x,  \
                             tk_z_y,  \
                             tk_z_z,  \
                             tk_zz_x, \
                             tk_zz_y, \
                             tk_zz_z, \
                             ts_0_x,  \
                             ts_0_y,  \
                             ts_0_z,  \
                             ts_zz_x, \
                             ts_zz_y, \
                             ts_zz_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_zz_x[i] = -2.0 * ts_0_x[i] * fbe_0 * fz_0 + tk_0_x[i] * fe_0 + tk_z_x[i] * pa_z[i] + 2.0 * ts_zz_x[i] * fz_0;

        tk_zz_y[i] = -2.0 * ts_0_y[i] * fbe_0 * fz_0 + tk_0_y[i] * fe_0 + tk_z_y[i] * pa_z[i] + 2.0 * ts_zz_y[i] * fz_0;

        tk_zz_z[i] = -2.0 * ts_0_z[i] * fbe_0 * fz_0 + tk_0_z[i] * fe_0 + tk_z_0[i] * fe_0 + tk_z_z[i] * pa_z[i] + 2.0 * ts_zz_z[i] * fz_0;
    }
}

}  // namespace kinrec
