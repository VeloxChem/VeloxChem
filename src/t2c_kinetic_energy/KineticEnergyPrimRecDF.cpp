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

#include "KineticEnergyPrimRecDF.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_df(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_df,
                            const size_t              idx_ovl_sf,
                            const size_t              idx_kin_sf,
                            const size_t              idx_kin_pd,
                            const size_t              idx_kin_pf,
                            const size_t              idx_ovl_df,
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

    // Set up components of auxiliary buffer : SF

    auto ts_0_xxx = pbuffer.data(idx_ovl_sf);

    auto ts_0_xxy = pbuffer.data(idx_ovl_sf + 1);

    auto ts_0_xxz = pbuffer.data(idx_ovl_sf + 2);

    auto ts_0_xyy = pbuffer.data(idx_ovl_sf + 3);

    auto ts_0_xyz = pbuffer.data(idx_ovl_sf + 4);

    auto ts_0_xzz = pbuffer.data(idx_ovl_sf + 5);

    auto ts_0_yyy = pbuffer.data(idx_ovl_sf + 6);

    auto ts_0_yyz = pbuffer.data(idx_ovl_sf + 7);

    auto ts_0_yzz = pbuffer.data(idx_ovl_sf + 8);

    auto ts_0_zzz = pbuffer.data(idx_ovl_sf + 9);

    // Set up components of auxiliary buffer : SF

    auto tk_0_xxx = pbuffer.data(idx_kin_sf);

    auto tk_0_xxy = pbuffer.data(idx_kin_sf + 1);

    auto tk_0_xxz = pbuffer.data(idx_kin_sf + 2);

    auto tk_0_xyy = pbuffer.data(idx_kin_sf + 3);

    auto tk_0_xyz = pbuffer.data(idx_kin_sf + 4);

    auto tk_0_xzz = pbuffer.data(idx_kin_sf + 5);

    auto tk_0_yyy = pbuffer.data(idx_kin_sf + 6);

    auto tk_0_yyz = pbuffer.data(idx_kin_sf + 7);

    auto tk_0_yzz = pbuffer.data(idx_kin_sf + 8);

    auto tk_0_zzz = pbuffer.data(idx_kin_sf + 9);

    // Set up components of auxiliary buffer : PD

    auto tk_x_xx = pbuffer.data(idx_kin_pd);

    auto tk_x_xy = pbuffer.data(idx_kin_pd + 1);

    auto tk_x_xz = pbuffer.data(idx_kin_pd + 2);

    auto tk_x_yy = pbuffer.data(idx_kin_pd + 3);

    auto tk_x_yz = pbuffer.data(idx_kin_pd + 4);

    auto tk_x_zz = pbuffer.data(idx_kin_pd + 5);

    auto tk_y_xx = pbuffer.data(idx_kin_pd + 6);

    auto tk_y_xy = pbuffer.data(idx_kin_pd + 7);

    auto tk_y_xz = pbuffer.data(idx_kin_pd + 8);

    auto tk_y_yy = pbuffer.data(idx_kin_pd + 9);

    auto tk_y_yz = pbuffer.data(idx_kin_pd + 10);

    auto tk_y_zz = pbuffer.data(idx_kin_pd + 11);

    auto tk_z_xx = pbuffer.data(idx_kin_pd + 12);

    auto tk_z_xy = pbuffer.data(idx_kin_pd + 13);

    auto tk_z_xz = pbuffer.data(idx_kin_pd + 14);

    auto tk_z_yy = pbuffer.data(idx_kin_pd + 15);

    auto tk_z_yz = pbuffer.data(idx_kin_pd + 16);

    auto tk_z_zz = pbuffer.data(idx_kin_pd + 17);

    // Set up components of auxiliary buffer : PF

    auto tk_x_xxx = pbuffer.data(idx_kin_pf);

    auto tk_x_xxy = pbuffer.data(idx_kin_pf + 1);

    auto tk_x_xxz = pbuffer.data(idx_kin_pf + 2);

    auto tk_x_xyy = pbuffer.data(idx_kin_pf + 3);

    auto tk_x_xyz = pbuffer.data(idx_kin_pf + 4);

    auto tk_x_xzz = pbuffer.data(idx_kin_pf + 5);

    auto tk_x_yyy = pbuffer.data(idx_kin_pf + 6);

    auto tk_x_yyz = pbuffer.data(idx_kin_pf + 7);

    auto tk_x_yzz = pbuffer.data(idx_kin_pf + 8);

    auto tk_x_zzz = pbuffer.data(idx_kin_pf + 9);

    auto tk_y_xxx = pbuffer.data(idx_kin_pf + 10);

    auto tk_y_xxy = pbuffer.data(idx_kin_pf + 11);

    auto tk_y_xxz = pbuffer.data(idx_kin_pf + 12);

    auto tk_y_xyy = pbuffer.data(idx_kin_pf + 13);

    auto tk_y_xyz = pbuffer.data(idx_kin_pf + 14);

    auto tk_y_xzz = pbuffer.data(idx_kin_pf + 15);

    auto tk_y_yyy = pbuffer.data(idx_kin_pf + 16);

    auto tk_y_yyz = pbuffer.data(idx_kin_pf + 17);

    auto tk_y_yzz = pbuffer.data(idx_kin_pf + 18);

    auto tk_y_zzz = pbuffer.data(idx_kin_pf + 19);

    auto tk_z_xxx = pbuffer.data(idx_kin_pf + 20);

    auto tk_z_xxy = pbuffer.data(idx_kin_pf + 21);

    auto tk_z_xxz = pbuffer.data(idx_kin_pf + 22);

    auto tk_z_xyy = pbuffer.data(idx_kin_pf + 23);

    auto tk_z_xyz = pbuffer.data(idx_kin_pf + 24);

    auto tk_z_xzz = pbuffer.data(idx_kin_pf + 25);

    auto tk_z_yyy = pbuffer.data(idx_kin_pf + 26);

    auto tk_z_yyz = pbuffer.data(idx_kin_pf + 27);

    auto tk_z_yzz = pbuffer.data(idx_kin_pf + 28);

    auto tk_z_zzz = pbuffer.data(idx_kin_pf + 29);

    // Set up components of auxiliary buffer : DF

    auto ts_xx_xxx = pbuffer.data(idx_ovl_df);

    auto ts_xx_xxy = pbuffer.data(idx_ovl_df + 1);

    auto ts_xx_xxz = pbuffer.data(idx_ovl_df + 2);

    auto ts_xx_xyy = pbuffer.data(idx_ovl_df + 3);

    auto ts_xx_xyz = pbuffer.data(idx_ovl_df + 4);

    auto ts_xx_xzz = pbuffer.data(idx_ovl_df + 5);

    auto ts_xx_yyy = pbuffer.data(idx_ovl_df + 6);

    auto ts_xx_yyz = pbuffer.data(idx_ovl_df + 7);

    auto ts_xx_yzz = pbuffer.data(idx_ovl_df + 8);

    auto ts_xx_zzz = pbuffer.data(idx_ovl_df + 9);

    auto ts_xy_xxx = pbuffer.data(idx_ovl_df + 10);

    auto ts_xy_xxy = pbuffer.data(idx_ovl_df + 11);

    auto ts_xy_xxz = pbuffer.data(idx_ovl_df + 12);

    auto ts_xy_xyy = pbuffer.data(idx_ovl_df + 13);

    auto ts_xy_xyz = pbuffer.data(idx_ovl_df + 14);

    auto ts_xy_xzz = pbuffer.data(idx_ovl_df + 15);

    auto ts_xy_yyy = pbuffer.data(idx_ovl_df + 16);

    auto ts_xy_yyz = pbuffer.data(idx_ovl_df + 17);

    auto ts_xy_yzz = pbuffer.data(idx_ovl_df + 18);

    auto ts_xy_zzz = pbuffer.data(idx_ovl_df + 19);

    auto ts_xz_xxx = pbuffer.data(idx_ovl_df + 20);

    auto ts_xz_xxy = pbuffer.data(idx_ovl_df + 21);

    auto ts_xz_xxz = pbuffer.data(idx_ovl_df + 22);

    auto ts_xz_xyy = pbuffer.data(idx_ovl_df + 23);

    auto ts_xz_xyz = pbuffer.data(idx_ovl_df + 24);

    auto ts_xz_xzz = pbuffer.data(idx_ovl_df + 25);

    auto ts_xz_yyy = pbuffer.data(idx_ovl_df + 26);

    auto ts_xz_yyz = pbuffer.data(idx_ovl_df + 27);

    auto ts_xz_yzz = pbuffer.data(idx_ovl_df + 28);

    auto ts_xz_zzz = pbuffer.data(idx_ovl_df + 29);

    auto ts_yy_xxx = pbuffer.data(idx_ovl_df + 30);

    auto ts_yy_xxy = pbuffer.data(idx_ovl_df + 31);

    auto ts_yy_xxz = pbuffer.data(idx_ovl_df + 32);

    auto ts_yy_xyy = pbuffer.data(idx_ovl_df + 33);

    auto ts_yy_xyz = pbuffer.data(idx_ovl_df + 34);

    auto ts_yy_xzz = pbuffer.data(idx_ovl_df + 35);

    auto ts_yy_yyy = pbuffer.data(idx_ovl_df + 36);

    auto ts_yy_yyz = pbuffer.data(idx_ovl_df + 37);

    auto ts_yy_yzz = pbuffer.data(idx_ovl_df + 38);

    auto ts_yy_zzz = pbuffer.data(idx_ovl_df + 39);

    auto ts_yz_xxx = pbuffer.data(idx_ovl_df + 40);

    auto ts_yz_xxy = pbuffer.data(idx_ovl_df + 41);

    auto ts_yz_xxz = pbuffer.data(idx_ovl_df + 42);

    auto ts_yz_xyy = pbuffer.data(idx_ovl_df + 43);

    auto ts_yz_xyz = pbuffer.data(idx_ovl_df + 44);

    auto ts_yz_xzz = pbuffer.data(idx_ovl_df + 45);

    auto ts_yz_yyy = pbuffer.data(idx_ovl_df + 46);

    auto ts_yz_yyz = pbuffer.data(idx_ovl_df + 47);

    auto ts_yz_yzz = pbuffer.data(idx_ovl_df + 48);

    auto ts_yz_zzz = pbuffer.data(idx_ovl_df + 49);

    auto ts_zz_xxx = pbuffer.data(idx_ovl_df + 50);

    auto ts_zz_xxy = pbuffer.data(idx_ovl_df + 51);

    auto ts_zz_xxz = pbuffer.data(idx_ovl_df + 52);

    auto ts_zz_xyy = pbuffer.data(idx_ovl_df + 53);

    auto ts_zz_xyz = pbuffer.data(idx_ovl_df + 54);

    auto ts_zz_xzz = pbuffer.data(idx_ovl_df + 55);

    auto ts_zz_yyy = pbuffer.data(idx_ovl_df + 56);

    auto ts_zz_yyz = pbuffer.data(idx_ovl_df + 57);

    auto ts_zz_yzz = pbuffer.data(idx_ovl_df + 58);

    auto ts_zz_zzz = pbuffer.data(idx_ovl_df + 59);

    // Set up 0-10 components of targeted buffer : DF

    auto tk_xx_xxx = pbuffer.data(idx_kin_df);

    auto tk_xx_xxy = pbuffer.data(idx_kin_df + 1);

    auto tk_xx_xxz = pbuffer.data(idx_kin_df + 2);

    auto tk_xx_xyy = pbuffer.data(idx_kin_df + 3);

    auto tk_xx_xyz = pbuffer.data(idx_kin_df + 4);

    auto tk_xx_xzz = pbuffer.data(idx_kin_df + 5);

    auto tk_xx_yyy = pbuffer.data(idx_kin_df + 6);

    auto tk_xx_yyz = pbuffer.data(idx_kin_df + 7);

    auto tk_xx_yzz = pbuffer.data(idx_kin_df + 8);

    auto tk_xx_zzz = pbuffer.data(idx_kin_df + 9);

#pragma omp simd aligned(pa_x,          \
                             tk_0_xxx,  \
                             tk_0_xxy,  \
                             tk_0_xxz,  \
                             tk_0_xyy,  \
                             tk_0_xyz,  \
                             tk_0_xzz,  \
                             tk_0_yyy,  \
                             tk_0_yyz,  \
                             tk_0_yzz,  \
                             tk_0_zzz,  \
                             tk_x_xx,   \
                             tk_x_xxx,  \
                             tk_x_xxy,  \
                             tk_x_xxz,  \
                             tk_x_xy,   \
                             tk_x_xyy,  \
                             tk_x_xyz,  \
                             tk_x_xz,   \
                             tk_x_xzz,  \
                             tk_x_yy,   \
                             tk_x_yyy,  \
                             tk_x_yyz,  \
                             tk_x_yz,   \
                             tk_x_yzz,  \
                             tk_x_zz,   \
                             tk_x_zzz,  \
                             tk_xx_xxx, \
                             tk_xx_xxy, \
                             tk_xx_xxz, \
                             tk_xx_xyy, \
                             tk_xx_xyz, \
                             tk_xx_xzz, \
                             tk_xx_yyy, \
                             tk_xx_yyz, \
                             tk_xx_yzz, \
                             tk_xx_zzz, \
                             ts_0_xxx,  \
                             ts_0_xxy,  \
                             ts_0_xxz,  \
                             ts_0_xyy,  \
                             ts_0_xyz,  \
                             ts_0_xzz,  \
                             ts_0_yyy,  \
                             ts_0_yyz,  \
                             ts_0_yzz,  \
                             ts_0_zzz,  \
                             ts_xx_xxx, \
                             ts_xx_xxy, \
                             ts_xx_xxz, \
                             ts_xx_xyy, \
                             ts_xx_xyz, \
                             ts_xx_xzz, \
                             ts_xx_yyy, \
                             ts_xx_yyz, \
                             ts_xx_yzz, \
                             ts_xx_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xx_xxx[i] =
            -2.0 * ts_0_xxx[i] * fbe_0 * fz_0 + tk_0_xxx[i] * fe_0 + 3.0 * tk_x_xx[i] * fe_0 + tk_x_xxx[i] * pa_x[i] + 2.0 * ts_xx_xxx[i] * fz_0;

        tk_xx_xxy[i] =
            -2.0 * ts_0_xxy[i] * fbe_0 * fz_0 + tk_0_xxy[i] * fe_0 + 2.0 * tk_x_xy[i] * fe_0 + tk_x_xxy[i] * pa_x[i] + 2.0 * ts_xx_xxy[i] * fz_0;

        tk_xx_xxz[i] =
            -2.0 * ts_0_xxz[i] * fbe_0 * fz_0 + tk_0_xxz[i] * fe_0 + 2.0 * tk_x_xz[i] * fe_0 + tk_x_xxz[i] * pa_x[i] + 2.0 * ts_xx_xxz[i] * fz_0;

        tk_xx_xyy[i] = -2.0 * ts_0_xyy[i] * fbe_0 * fz_0 + tk_0_xyy[i] * fe_0 + tk_x_yy[i] * fe_0 + tk_x_xyy[i] * pa_x[i] + 2.0 * ts_xx_xyy[i] * fz_0;

        tk_xx_xyz[i] = -2.0 * ts_0_xyz[i] * fbe_0 * fz_0 + tk_0_xyz[i] * fe_0 + tk_x_yz[i] * fe_0 + tk_x_xyz[i] * pa_x[i] + 2.0 * ts_xx_xyz[i] * fz_0;

        tk_xx_xzz[i] = -2.0 * ts_0_xzz[i] * fbe_0 * fz_0 + tk_0_xzz[i] * fe_0 + tk_x_zz[i] * fe_0 + tk_x_xzz[i] * pa_x[i] + 2.0 * ts_xx_xzz[i] * fz_0;

        tk_xx_yyy[i] = -2.0 * ts_0_yyy[i] * fbe_0 * fz_0 + tk_0_yyy[i] * fe_0 + tk_x_yyy[i] * pa_x[i] + 2.0 * ts_xx_yyy[i] * fz_0;

        tk_xx_yyz[i] = -2.0 * ts_0_yyz[i] * fbe_0 * fz_0 + tk_0_yyz[i] * fe_0 + tk_x_yyz[i] * pa_x[i] + 2.0 * ts_xx_yyz[i] * fz_0;

        tk_xx_yzz[i] = -2.0 * ts_0_yzz[i] * fbe_0 * fz_0 + tk_0_yzz[i] * fe_0 + tk_x_yzz[i] * pa_x[i] + 2.0 * ts_xx_yzz[i] * fz_0;

        tk_xx_zzz[i] = -2.0 * ts_0_zzz[i] * fbe_0 * fz_0 + tk_0_zzz[i] * fe_0 + tk_x_zzz[i] * pa_x[i] + 2.0 * ts_xx_zzz[i] * fz_0;
    }

    // Set up 10-20 components of targeted buffer : DF

    auto tk_xy_xxx = pbuffer.data(idx_kin_df + 10);

    auto tk_xy_xxy = pbuffer.data(idx_kin_df + 11);

    auto tk_xy_xxz = pbuffer.data(idx_kin_df + 12);

    auto tk_xy_xyy = pbuffer.data(idx_kin_df + 13);

    auto tk_xy_xyz = pbuffer.data(idx_kin_df + 14);

    auto tk_xy_xzz = pbuffer.data(idx_kin_df + 15);

    auto tk_xy_yyy = pbuffer.data(idx_kin_df + 16);

    auto tk_xy_yyz = pbuffer.data(idx_kin_df + 17);

    auto tk_xy_yzz = pbuffer.data(idx_kin_df + 18);

    auto tk_xy_zzz = pbuffer.data(idx_kin_df + 19);

#pragma omp simd aligned(pa_x,          \
                             pa_y,      \
                             tk_x_xxx,  \
                             tk_x_xxz,  \
                             tk_x_xzz,  \
                             tk_xy_xxx, \
                             tk_xy_xxy, \
                             tk_xy_xxz, \
                             tk_xy_xyy, \
                             tk_xy_xyz, \
                             tk_xy_xzz, \
                             tk_xy_yyy, \
                             tk_xy_yyz, \
                             tk_xy_yzz, \
                             tk_xy_zzz, \
                             tk_y_xxy,  \
                             tk_y_xy,   \
                             tk_y_xyy,  \
                             tk_y_xyz,  \
                             tk_y_yy,   \
                             tk_y_yyy,  \
                             tk_y_yyz,  \
                             tk_y_yz,   \
                             tk_y_yzz,  \
                             tk_y_zzz,  \
                             ts_xy_xxx, \
                             ts_xy_xxy, \
                             ts_xy_xxz, \
                             ts_xy_xyy, \
                             ts_xy_xyz, \
                             ts_xy_xzz, \
                             ts_xy_yyy, \
                             ts_xy_yyz, \
                             ts_xy_yzz, \
                             ts_xy_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xy_xxx[i] = tk_x_xxx[i] * pa_y[i] + 2.0 * ts_xy_xxx[i] * fz_0;

        tk_xy_xxy[i] = 2.0 * tk_y_xy[i] * fe_0 + tk_y_xxy[i] * pa_x[i] + 2.0 * ts_xy_xxy[i] * fz_0;

        tk_xy_xxz[i] = tk_x_xxz[i] * pa_y[i] + 2.0 * ts_xy_xxz[i] * fz_0;

        tk_xy_xyy[i] = tk_y_yy[i] * fe_0 + tk_y_xyy[i] * pa_x[i] + 2.0 * ts_xy_xyy[i] * fz_0;

        tk_xy_xyz[i] = tk_y_yz[i] * fe_0 + tk_y_xyz[i] * pa_x[i] + 2.0 * ts_xy_xyz[i] * fz_0;

        tk_xy_xzz[i] = tk_x_xzz[i] * pa_y[i] + 2.0 * ts_xy_xzz[i] * fz_0;

        tk_xy_yyy[i] = tk_y_yyy[i] * pa_x[i] + 2.0 * ts_xy_yyy[i] * fz_0;

        tk_xy_yyz[i] = tk_y_yyz[i] * pa_x[i] + 2.0 * ts_xy_yyz[i] * fz_0;

        tk_xy_yzz[i] = tk_y_yzz[i] * pa_x[i] + 2.0 * ts_xy_yzz[i] * fz_0;

        tk_xy_zzz[i] = tk_y_zzz[i] * pa_x[i] + 2.0 * ts_xy_zzz[i] * fz_0;
    }

    // Set up 20-30 components of targeted buffer : DF

    auto tk_xz_xxx = pbuffer.data(idx_kin_df + 20);

    auto tk_xz_xxy = pbuffer.data(idx_kin_df + 21);

    auto tk_xz_xxz = pbuffer.data(idx_kin_df + 22);

    auto tk_xz_xyy = pbuffer.data(idx_kin_df + 23);

    auto tk_xz_xyz = pbuffer.data(idx_kin_df + 24);

    auto tk_xz_xzz = pbuffer.data(idx_kin_df + 25);

    auto tk_xz_yyy = pbuffer.data(idx_kin_df + 26);

    auto tk_xz_yyz = pbuffer.data(idx_kin_df + 27);

    auto tk_xz_yzz = pbuffer.data(idx_kin_df + 28);

    auto tk_xz_zzz = pbuffer.data(idx_kin_df + 29);

#pragma omp simd aligned(pa_x,          \
                             pa_z,      \
                             tk_x_xxx,  \
                             tk_x_xxy,  \
                             tk_x_xyy,  \
                             tk_xz_xxx, \
                             tk_xz_xxy, \
                             tk_xz_xxz, \
                             tk_xz_xyy, \
                             tk_xz_xyz, \
                             tk_xz_xzz, \
                             tk_xz_yyy, \
                             tk_xz_yyz, \
                             tk_xz_yzz, \
                             tk_xz_zzz, \
                             tk_z_xxz,  \
                             tk_z_xyz,  \
                             tk_z_xz,   \
                             tk_z_xzz,  \
                             tk_z_yyy,  \
                             tk_z_yyz,  \
                             tk_z_yz,   \
                             tk_z_yzz,  \
                             tk_z_zz,   \
                             tk_z_zzz,  \
                             ts_xz_xxx, \
                             ts_xz_xxy, \
                             ts_xz_xxz, \
                             ts_xz_xyy, \
                             ts_xz_xyz, \
                             ts_xz_xzz, \
                             ts_xz_yyy, \
                             ts_xz_yyz, \
                             ts_xz_yzz, \
                             ts_xz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xz_xxx[i] = tk_x_xxx[i] * pa_z[i] + 2.0 * ts_xz_xxx[i] * fz_0;

        tk_xz_xxy[i] = tk_x_xxy[i] * pa_z[i] + 2.0 * ts_xz_xxy[i] * fz_0;

        tk_xz_xxz[i] = 2.0 * tk_z_xz[i] * fe_0 + tk_z_xxz[i] * pa_x[i] + 2.0 * ts_xz_xxz[i] * fz_0;

        tk_xz_xyy[i] = tk_x_xyy[i] * pa_z[i] + 2.0 * ts_xz_xyy[i] * fz_0;

        tk_xz_xyz[i] = tk_z_yz[i] * fe_0 + tk_z_xyz[i] * pa_x[i] + 2.0 * ts_xz_xyz[i] * fz_0;

        tk_xz_xzz[i] = tk_z_zz[i] * fe_0 + tk_z_xzz[i] * pa_x[i] + 2.0 * ts_xz_xzz[i] * fz_0;

        tk_xz_yyy[i] = tk_z_yyy[i] * pa_x[i] + 2.0 * ts_xz_yyy[i] * fz_0;

        tk_xz_yyz[i] = tk_z_yyz[i] * pa_x[i] + 2.0 * ts_xz_yyz[i] * fz_0;

        tk_xz_yzz[i] = tk_z_yzz[i] * pa_x[i] + 2.0 * ts_xz_yzz[i] * fz_0;

        tk_xz_zzz[i] = tk_z_zzz[i] * pa_x[i] + 2.0 * ts_xz_zzz[i] * fz_0;
    }

    // Set up 30-40 components of targeted buffer : DF

    auto tk_yy_xxx = pbuffer.data(idx_kin_df + 30);

    auto tk_yy_xxy = pbuffer.data(idx_kin_df + 31);

    auto tk_yy_xxz = pbuffer.data(idx_kin_df + 32);

    auto tk_yy_xyy = pbuffer.data(idx_kin_df + 33);

    auto tk_yy_xyz = pbuffer.data(idx_kin_df + 34);

    auto tk_yy_xzz = pbuffer.data(idx_kin_df + 35);

    auto tk_yy_yyy = pbuffer.data(idx_kin_df + 36);

    auto tk_yy_yyz = pbuffer.data(idx_kin_df + 37);

    auto tk_yy_yzz = pbuffer.data(idx_kin_df + 38);

    auto tk_yy_zzz = pbuffer.data(idx_kin_df + 39);

#pragma omp simd aligned(pa_y,          \
                             tk_0_xxx,  \
                             tk_0_xxy,  \
                             tk_0_xxz,  \
                             tk_0_xyy,  \
                             tk_0_xyz,  \
                             tk_0_xzz,  \
                             tk_0_yyy,  \
                             tk_0_yyz,  \
                             tk_0_yzz,  \
                             tk_0_zzz,  \
                             tk_y_xx,   \
                             tk_y_xxx,  \
                             tk_y_xxy,  \
                             tk_y_xxz,  \
                             tk_y_xy,   \
                             tk_y_xyy,  \
                             tk_y_xyz,  \
                             tk_y_xz,   \
                             tk_y_xzz,  \
                             tk_y_yy,   \
                             tk_y_yyy,  \
                             tk_y_yyz,  \
                             tk_y_yz,   \
                             tk_y_yzz,  \
                             tk_y_zz,   \
                             tk_y_zzz,  \
                             tk_yy_xxx, \
                             tk_yy_xxy, \
                             tk_yy_xxz, \
                             tk_yy_xyy, \
                             tk_yy_xyz, \
                             tk_yy_xzz, \
                             tk_yy_yyy, \
                             tk_yy_yyz, \
                             tk_yy_yzz, \
                             tk_yy_zzz, \
                             ts_0_xxx,  \
                             ts_0_xxy,  \
                             ts_0_xxz,  \
                             ts_0_xyy,  \
                             ts_0_xyz,  \
                             ts_0_xzz,  \
                             ts_0_yyy,  \
                             ts_0_yyz,  \
                             ts_0_yzz,  \
                             ts_0_zzz,  \
                             ts_yy_xxx, \
                             ts_yy_xxy, \
                             ts_yy_xxz, \
                             ts_yy_xyy, \
                             ts_yy_xyz, \
                             ts_yy_xzz, \
                             ts_yy_yyy, \
                             ts_yy_yyz, \
                             ts_yy_yzz, \
                             ts_yy_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yy_xxx[i] = -2.0 * ts_0_xxx[i] * fbe_0 * fz_0 + tk_0_xxx[i] * fe_0 + tk_y_xxx[i] * pa_y[i] + 2.0 * ts_yy_xxx[i] * fz_0;

        tk_yy_xxy[i] = -2.0 * ts_0_xxy[i] * fbe_0 * fz_0 + tk_0_xxy[i] * fe_0 + tk_y_xx[i] * fe_0 + tk_y_xxy[i] * pa_y[i] + 2.0 * ts_yy_xxy[i] * fz_0;

        tk_yy_xxz[i] = -2.0 * ts_0_xxz[i] * fbe_0 * fz_0 + tk_0_xxz[i] * fe_0 + tk_y_xxz[i] * pa_y[i] + 2.0 * ts_yy_xxz[i] * fz_0;

        tk_yy_xyy[i] =
            -2.0 * ts_0_xyy[i] * fbe_0 * fz_0 + tk_0_xyy[i] * fe_0 + 2.0 * tk_y_xy[i] * fe_0 + tk_y_xyy[i] * pa_y[i] + 2.0 * ts_yy_xyy[i] * fz_0;

        tk_yy_xyz[i] = -2.0 * ts_0_xyz[i] * fbe_0 * fz_0 + tk_0_xyz[i] * fe_0 + tk_y_xz[i] * fe_0 + tk_y_xyz[i] * pa_y[i] + 2.0 * ts_yy_xyz[i] * fz_0;

        tk_yy_xzz[i] = -2.0 * ts_0_xzz[i] * fbe_0 * fz_0 + tk_0_xzz[i] * fe_0 + tk_y_xzz[i] * pa_y[i] + 2.0 * ts_yy_xzz[i] * fz_0;

        tk_yy_yyy[i] =
            -2.0 * ts_0_yyy[i] * fbe_0 * fz_0 + tk_0_yyy[i] * fe_0 + 3.0 * tk_y_yy[i] * fe_0 + tk_y_yyy[i] * pa_y[i] + 2.0 * ts_yy_yyy[i] * fz_0;

        tk_yy_yyz[i] =
            -2.0 * ts_0_yyz[i] * fbe_0 * fz_0 + tk_0_yyz[i] * fe_0 + 2.0 * tk_y_yz[i] * fe_0 + tk_y_yyz[i] * pa_y[i] + 2.0 * ts_yy_yyz[i] * fz_0;

        tk_yy_yzz[i] = -2.0 * ts_0_yzz[i] * fbe_0 * fz_0 + tk_0_yzz[i] * fe_0 + tk_y_zz[i] * fe_0 + tk_y_yzz[i] * pa_y[i] + 2.0 * ts_yy_yzz[i] * fz_0;

        tk_yy_zzz[i] = -2.0 * ts_0_zzz[i] * fbe_0 * fz_0 + tk_0_zzz[i] * fe_0 + tk_y_zzz[i] * pa_y[i] + 2.0 * ts_yy_zzz[i] * fz_0;
    }

    // Set up 40-50 components of targeted buffer : DF

    auto tk_yz_xxx = pbuffer.data(idx_kin_df + 40);

    auto tk_yz_xxy = pbuffer.data(idx_kin_df + 41);

    auto tk_yz_xxz = pbuffer.data(idx_kin_df + 42);

    auto tk_yz_xyy = pbuffer.data(idx_kin_df + 43);

    auto tk_yz_xyz = pbuffer.data(idx_kin_df + 44);

    auto tk_yz_xzz = pbuffer.data(idx_kin_df + 45);

    auto tk_yz_yyy = pbuffer.data(idx_kin_df + 46);

    auto tk_yz_yyz = pbuffer.data(idx_kin_df + 47);

    auto tk_yz_yzz = pbuffer.data(idx_kin_df + 48);

    auto tk_yz_zzz = pbuffer.data(idx_kin_df + 49);

#pragma omp simd aligned(pa_y,          \
                             pa_z,      \
                             tk_y_xxy,  \
                             tk_y_xyy,  \
                             tk_y_yyy,  \
                             tk_yz_xxx, \
                             tk_yz_xxy, \
                             tk_yz_xxz, \
                             tk_yz_xyy, \
                             tk_yz_xyz, \
                             tk_yz_xzz, \
                             tk_yz_yyy, \
                             tk_yz_yyz, \
                             tk_yz_yzz, \
                             tk_yz_zzz, \
                             tk_z_xxx,  \
                             tk_z_xxz,  \
                             tk_z_xyz,  \
                             tk_z_xz,   \
                             tk_z_xzz,  \
                             tk_z_yyz,  \
                             tk_z_yz,   \
                             tk_z_yzz,  \
                             tk_z_zz,   \
                             tk_z_zzz,  \
                             ts_yz_xxx, \
                             ts_yz_xxy, \
                             ts_yz_xxz, \
                             ts_yz_xyy, \
                             ts_yz_xyz, \
                             ts_yz_xzz, \
                             ts_yz_yyy, \
                             ts_yz_yyz, \
                             ts_yz_yzz, \
                             ts_yz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yz_xxx[i] = tk_z_xxx[i] * pa_y[i] + 2.0 * ts_yz_xxx[i] * fz_0;

        tk_yz_xxy[i] = tk_y_xxy[i] * pa_z[i] + 2.0 * ts_yz_xxy[i] * fz_0;

        tk_yz_xxz[i] = tk_z_xxz[i] * pa_y[i] + 2.0 * ts_yz_xxz[i] * fz_0;

        tk_yz_xyy[i] = tk_y_xyy[i] * pa_z[i] + 2.0 * ts_yz_xyy[i] * fz_0;

        tk_yz_xyz[i] = tk_z_xz[i] * fe_0 + tk_z_xyz[i] * pa_y[i] + 2.0 * ts_yz_xyz[i] * fz_0;

        tk_yz_xzz[i] = tk_z_xzz[i] * pa_y[i] + 2.0 * ts_yz_xzz[i] * fz_0;

        tk_yz_yyy[i] = tk_y_yyy[i] * pa_z[i] + 2.0 * ts_yz_yyy[i] * fz_0;

        tk_yz_yyz[i] = 2.0 * tk_z_yz[i] * fe_0 + tk_z_yyz[i] * pa_y[i] + 2.0 * ts_yz_yyz[i] * fz_0;

        tk_yz_yzz[i] = tk_z_zz[i] * fe_0 + tk_z_yzz[i] * pa_y[i] + 2.0 * ts_yz_yzz[i] * fz_0;

        tk_yz_zzz[i] = tk_z_zzz[i] * pa_y[i] + 2.0 * ts_yz_zzz[i] * fz_0;
    }

    // Set up 50-60 components of targeted buffer : DF

    auto tk_zz_xxx = pbuffer.data(idx_kin_df + 50);

    auto tk_zz_xxy = pbuffer.data(idx_kin_df + 51);

    auto tk_zz_xxz = pbuffer.data(idx_kin_df + 52);

    auto tk_zz_xyy = pbuffer.data(idx_kin_df + 53);

    auto tk_zz_xyz = pbuffer.data(idx_kin_df + 54);

    auto tk_zz_xzz = pbuffer.data(idx_kin_df + 55);

    auto tk_zz_yyy = pbuffer.data(idx_kin_df + 56);

    auto tk_zz_yyz = pbuffer.data(idx_kin_df + 57);

    auto tk_zz_yzz = pbuffer.data(idx_kin_df + 58);

    auto tk_zz_zzz = pbuffer.data(idx_kin_df + 59);

#pragma omp simd aligned(pa_z,          \
                             tk_0_xxx,  \
                             tk_0_xxy,  \
                             tk_0_xxz,  \
                             tk_0_xyy,  \
                             tk_0_xyz,  \
                             tk_0_xzz,  \
                             tk_0_yyy,  \
                             tk_0_yyz,  \
                             tk_0_yzz,  \
                             tk_0_zzz,  \
                             tk_z_xx,   \
                             tk_z_xxx,  \
                             tk_z_xxy,  \
                             tk_z_xxz,  \
                             tk_z_xy,   \
                             tk_z_xyy,  \
                             tk_z_xyz,  \
                             tk_z_xz,   \
                             tk_z_xzz,  \
                             tk_z_yy,   \
                             tk_z_yyy,  \
                             tk_z_yyz,  \
                             tk_z_yz,   \
                             tk_z_yzz,  \
                             tk_z_zz,   \
                             tk_z_zzz,  \
                             tk_zz_xxx, \
                             tk_zz_xxy, \
                             tk_zz_xxz, \
                             tk_zz_xyy, \
                             tk_zz_xyz, \
                             tk_zz_xzz, \
                             tk_zz_yyy, \
                             tk_zz_yyz, \
                             tk_zz_yzz, \
                             tk_zz_zzz, \
                             ts_0_xxx,  \
                             ts_0_xxy,  \
                             ts_0_xxz,  \
                             ts_0_xyy,  \
                             ts_0_xyz,  \
                             ts_0_xzz,  \
                             ts_0_yyy,  \
                             ts_0_yyz,  \
                             ts_0_yzz,  \
                             ts_0_zzz,  \
                             ts_zz_xxx, \
                             ts_zz_xxy, \
                             ts_zz_xxz, \
                             ts_zz_xyy, \
                             ts_zz_xyz, \
                             ts_zz_xzz, \
                             ts_zz_yyy, \
                             ts_zz_yyz, \
                             ts_zz_yzz, \
                             ts_zz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_zz_xxx[i] = -2.0 * ts_0_xxx[i] * fbe_0 * fz_0 + tk_0_xxx[i] * fe_0 + tk_z_xxx[i] * pa_z[i] + 2.0 * ts_zz_xxx[i] * fz_0;

        tk_zz_xxy[i] = -2.0 * ts_0_xxy[i] * fbe_0 * fz_0 + tk_0_xxy[i] * fe_0 + tk_z_xxy[i] * pa_z[i] + 2.0 * ts_zz_xxy[i] * fz_0;

        tk_zz_xxz[i] = -2.0 * ts_0_xxz[i] * fbe_0 * fz_0 + tk_0_xxz[i] * fe_0 + tk_z_xx[i] * fe_0 + tk_z_xxz[i] * pa_z[i] + 2.0 * ts_zz_xxz[i] * fz_0;

        tk_zz_xyy[i] = -2.0 * ts_0_xyy[i] * fbe_0 * fz_0 + tk_0_xyy[i] * fe_0 + tk_z_xyy[i] * pa_z[i] + 2.0 * ts_zz_xyy[i] * fz_0;

        tk_zz_xyz[i] = -2.0 * ts_0_xyz[i] * fbe_0 * fz_0 + tk_0_xyz[i] * fe_0 + tk_z_xy[i] * fe_0 + tk_z_xyz[i] * pa_z[i] + 2.0 * ts_zz_xyz[i] * fz_0;

        tk_zz_xzz[i] =
            -2.0 * ts_0_xzz[i] * fbe_0 * fz_0 + tk_0_xzz[i] * fe_0 + 2.0 * tk_z_xz[i] * fe_0 + tk_z_xzz[i] * pa_z[i] + 2.0 * ts_zz_xzz[i] * fz_0;

        tk_zz_yyy[i] = -2.0 * ts_0_yyy[i] * fbe_0 * fz_0 + tk_0_yyy[i] * fe_0 + tk_z_yyy[i] * pa_z[i] + 2.0 * ts_zz_yyy[i] * fz_0;

        tk_zz_yyz[i] = -2.0 * ts_0_yyz[i] * fbe_0 * fz_0 + tk_0_yyz[i] * fe_0 + tk_z_yy[i] * fe_0 + tk_z_yyz[i] * pa_z[i] + 2.0 * ts_zz_yyz[i] * fz_0;

        tk_zz_yzz[i] =
            -2.0 * ts_0_yzz[i] * fbe_0 * fz_0 + tk_0_yzz[i] * fe_0 + 2.0 * tk_z_yz[i] * fe_0 + tk_z_yzz[i] * pa_z[i] + 2.0 * ts_zz_yzz[i] * fz_0;

        tk_zz_zzz[i] =
            -2.0 * ts_0_zzz[i] * fbe_0 * fz_0 + tk_0_zzz[i] * fe_0 + 3.0 * tk_z_zz[i] * fe_0 + tk_z_zzz[i] * pa_z[i] + 2.0 * ts_zz_zzz[i] * fz_0;
    }
}

}  // namespace kinrec
