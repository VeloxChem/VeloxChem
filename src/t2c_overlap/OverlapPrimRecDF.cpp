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

#include "OverlapPrimRecDF.hpp"

namespace ovlrec {  // ovlrec namespace

auto
comp_prim_overlap_df(CSimdArray<double>&       pbuffer,
                     const size_t              idx_ovl_df,
                     const size_t              idx_ovl_sf,
                     const size_t              idx_ovl_pd,
                     const size_t              idx_ovl_pf,
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

    // Set up components of auxiliary buffer : PD

    auto ts_x_xx = pbuffer.data(idx_ovl_pd);

    auto ts_x_xy = pbuffer.data(idx_ovl_pd + 1);

    auto ts_x_xz = pbuffer.data(idx_ovl_pd + 2);

    auto ts_x_yy = pbuffer.data(idx_ovl_pd + 3);

    auto ts_x_yz = pbuffer.data(idx_ovl_pd + 4);

    auto ts_x_zz = pbuffer.data(idx_ovl_pd + 5);

    auto ts_y_xx = pbuffer.data(idx_ovl_pd + 6);

    auto ts_y_xy = pbuffer.data(idx_ovl_pd + 7);

    auto ts_y_xz = pbuffer.data(idx_ovl_pd + 8);

    auto ts_y_yy = pbuffer.data(idx_ovl_pd + 9);

    auto ts_y_yz = pbuffer.data(idx_ovl_pd + 10);

    auto ts_y_zz = pbuffer.data(idx_ovl_pd + 11);

    auto ts_z_xx = pbuffer.data(idx_ovl_pd + 12);

    auto ts_z_xy = pbuffer.data(idx_ovl_pd + 13);

    auto ts_z_xz = pbuffer.data(idx_ovl_pd + 14);

    auto ts_z_yy = pbuffer.data(idx_ovl_pd + 15);

    auto ts_z_yz = pbuffer.data(idx_ovl_pd + 16);

    auto ts_z_zz = pbuffer.data(idx_ovl_pd + 17);

    // Set up components of auxiliary buffer : PF

    auto ts_x_xxx = pbuffer.data(idx_ovl_pf);

    auto ts_x_xxy = pbuffer.data(idx_ovl_pf + 1);

    auto ts_x_xxz = pbuffer.data(idx_ovl_pf + 2);

    auto ts_x_xyy = pbuffer.data(idx_ovl_pf + 3);

    auto ts_x_xyz = pbuffer.data(idx_ovl_pf + 4);

    auto ts_x_xzz = pbuffer.data(idx_ovl_pf + 5);

    auto ts_x_yyy = pbuffer.data(idx_ovl_pf + 6);

    auto ts_x_yyz = pbuffer.data(idx_ovl_pf + 7);

    auto ts_x_yzz = pbuffer.data(idx_ovl_pf + 8);

    auto ts_x_zzz = pbuffer.data(idx_ovl_pf + 9);

    auto ts_y_xxx = pbuffer.data(idx_ovl_pf + 10);

    auto ts_y_xxy = pbuffer.data(idx_ovl_pf + 11);

    auto ts_y_xxz = pbuffer.data(idx_ovl_pf + 12);

    auto ts_y_xyy = pbuffer.data(idx_ovl_pf + 13);

    auto ts_y_xyz = pbuffer.data(idx_ovl_pf + 14);

    auto ts_y_xzz = pbuffer.data(idx_ovl_pf + 15);

    auto ts_y_yyy = pbuffer.data(idx_ovl_pf + 16);

    auto ts_y_yyz = pbuffer.data(idx_ovl_pf + 17);

    auto ts_y_yzz = pbuffer.data(idx_ovl_pf + 18);

    auto ts_y_zzz = pbuffer.data(idx_ovl_pf + 19);

    auto ts_z_xxx = pbuffer.data(idx_ovl_pf + 20);

    auto ts_z_xxy = pbuffer.data(idx_ovl_pf + 21);

    auto ts_z_xxz = pbuffer.data(idx_ovl_pf + 22);

    auto ts_z_xyy = pbuffer.data(idx_ovl_pf + 23);

    auto ts_z_xyz = pbuffer.data(idx_ovl_pf + 24);

    auto ts_z_xzz = pbuffer.data(idx_ovl_pf + 25);

    auto ts_z_yyy = pbuffer.data(idx_ovl_pf + 26);

    auto ts_z_yyz = pbuffer.data(idx_ovl_pf + 27);

    auto ts_z_yzz = pbuffer.data(idx_ovl_pf + 28);

    auto ts_z_zzz = pbuffer.data(idx_ovl_pf + 29);

    // Set up 0-10 components of targeted buffer : DF

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

#pragma omp simd aligned(pa_x,          \
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
                             ts_x_xx,   \
                             ts_x_xxx,  \
                             ts_x_xxy,  \
                             ts_x_xxz,  \
                             ts_x_xy,   \
                             ts_x_xyy,  \
                             ts_x_xyz,  \
                             ts_x_xz,   \
                             ts_x_xzz,  \
                             ts_x_yy,   \
                             ts_x_yyy,  \
                             ts_x_yyz,  \
                             ts_x_yz,   \
                             ts_x_yzz,  \
                             ts_x_zz,   \
                             ts_x_zzz,  \
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

        ts_xx_xxx[i] = ts_0_xxx[i] * fe_0 + 3.0 * ts_x_xx[i] * fe_0 + ts_x_xxx[i] * pa_x[i];

        ts_xx_xxy[i] = ts_0_xxy[i] * fe_0 + 2.0 * ts_x_xy[i] * fe_0 + ts_x_xxy[i] * pa_x[i];

        ts_xx_xxz[i] = ts_0_xxz[i] * fe_0 + 2.0 * ts_x_xz[i] * fe_0 + ts_x_xxz[i] * pa_x[i];

        ts_xx_xyy[i] = ts_0_xyy[i] * fe_0 + ts_x_yy[i] * fe_0 + ts_x_xyy[i] * pa_x[i];

        ts_xx_xyz[i] = ts_0_xyz[i] * fe_0 + ts_x_yz[i] * fe_0 + ts_x_xyz[i] * pa_x[i];

        ts_xx_xzz[i] = ts_0_xzz[i] * fe_0 + ts_x_zz[i] * fe_0 + ts_x_xzz[i] * pa_x[i];

        ts_xx_yyy[i] = ts_0_yyy[i] * fe_0 + ts_x_yyy[i] * pa_x[i];

        ts_xx_yyz[i] = ts_0_yyz[i] * fe_0 + ts_x_yyz[i] * pa_x[i];

        ts_xx_yzz[i] = ts_0_yzz[i] * fe_0 + ts_x_yzz[i] * pa_x[i];

        ts_xx_zzz[i] = ts_0_zzz[i] * fe_0 + ts_x_zzz[i] * pa_x[i];
    }

    // Set up 10-20 components of targeted buffer : DF

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

#pragma omp simd aligned(pa_x,          \
                             pa_y,      \
                             ts_x_xxx,  \
                             ts_x_xxz,  \
                             ts_x_xzz,  \
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
                             ts_y_xxy,  \
                             ts_y_xy,   \
                             ts_y_xyy,  \
                             ts_y_xyz,  \
                             ts_y_yy,   \
                             ts_y_yyy,  \
                             ts_y_yyz,  \
                             ts_y_yz,   \
                             ts_y_yzz,  \
                             ts_y_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xy_xxx[i] = ts_x_xxx[i] * pa_y[i];

        ts_xy_xxy[i] = 2.0 * ts_y_xy[i] * fe_0 + ts_y_xxy[i] * pa_x[i];

        ts_xy_xxz[i] = ts_x_xxz[i] * pa_y[i];

        ts_xy_xyy[i] = ts_y_yy[i] * fe_0 + ts_y_xyy[i] * pa_x[i];

        ts_xy_xyz[i] = ts_y_yz[i] * fe_0 + ts_y_xyz[i] * pa_x[i];

        ts_xy_xzz[i] = ts_x_xzz[i] * pa_y[i];

        ts_xy_yyy[i] = ts_y_yyy[i] * pa_x[i];

        ts_xy_yyz[i] = ts_y_yyz[i] * pa_x[i];

        ts_xy_yzz[i] = ts_y_yzz[i] * pa_x[i];

        ts_xy_zzz[i] = ts_y_zzz[i] * pa_x[i];
    }

    // Set up 20-30 components of targeted buffer : DF

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

#pragma omp simd aligned(pa_x,          \
                             pa_z,      \
                             ts_x_xxx,  \
                             ts_x_xxy,  \
                             ts_x_xyy,  \
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
                             ts_z_xxz,  \
                             ts_z_xyz,  \
                             ts_z_xz,   \
                             ts_z_xzz,  \
                             ts_z_yyy,  \
                             ts_z_yyz,  \
                             ts_z_yz,   \
                             ts_z_yzz,  \
                             ts_z_zz,   \
                             ts_z_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xz_xxx[i] = ts_x_xxx[i] * pa_z[i];

        ts_xz_xxy[i] = ts_x_xxy[i] * pa_z[i];

        ts_xz_xxz[i] = 2.0 * ts_z_xz[i] * fe_0 + ts_z_xxz[i] * pa_x[i];

        ts_xz_xyy[i] = ts_x_xyy[i] * pa_z[i];

        ts_xz_xyz[i] = ts_z_yz[i] * fe_0 + ts_z_xyz[i] * pa_x[i];

        ts_xz_xzz[i] = ts_z_zz[i] * fe_0 + ts_z_xzz[i] * pa_x[i];

        ts_xz_yyy[i] = ts_z_yyy[i] * pa_x[i];

        ts_xz_yyz[i] = ts_z_yyz[i] * pa_x[i];

        ts_xz_yzz[i] = ts_z_yzz[i] * pa_x[i];

        ts_xz_zzz[i] = ts_z_zzz[i] * pa_x[i];
    }

    // Set up 30-40 components of targeted buffer : DF

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

#pragma omp simd aligned(pa_y,          \
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
                             ts_y_xx,   \
                             ts_y_xxx,  \
                             ts_y_xxy,  \
                             ts_y_xxz,  \
                             ts_y_xy,   \
                             ts_y_xyy,  \
                             ts_y_xyz,  \
                             ts_y_xz,   \
                             ts_y_xzz,  \
                             ts_y_yy,   \
                             ts_y_yyy,  \
                             ts_y_yyz,  \
                             ts_y_yz,   \
                             ts_y_yzz,  \
                             ts_y_zz,   \
                             ts_y_zzz,  \
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

        ts_yy_xxx[i] = ts_0_xxx[i] * fe_0 + ts_y_xxx[i] * pa_y[i];

        ts_yy_xxy[i] = ts_0_xxy[i] * fe_0 + ts_y_xx[i] * fe_0 + ts_y_xxy[i] * pa_y[i];

        ts_yy_xxz[i] = ts_0_xxz[i] * fe_0 + ts_y_xxz[i] * pa_y[i];

        ts_yy_xyy[i] = ts_0_xyy[i] * fe_0 + 2.0 * ts_y_xy[i] * fe_0 + ts_y_xyy[i] * pa_y[i];

        ts_yy_xyz[i] = ts_0_xyz[i] * fe_0 + ts_y_xz[i] * fe_0 + ts_y_xyz[i] * pa_y[i];

        ts_yy_xzz[i] = ts_0_xzz[i] * fe_0 + ts_y_xzz[i] * pa_y[i];

        ts_yy_yyy[i] = ts_0_yyy[i] * fe_0 + 3.0 * ts_y_yy[i] * fe_0 + ts_y_yyy[i] * pa_y[i];

        ts_yy_yyz[i] = ts_0_yyz[i] * fe_0 + 2.0 * ts_y_yz[i] * fe_0 + ts_y_yyz[i] * pa_y[i];

        ts_yy_yzz[i] = ts_0_yzz[i] * fe_0 + ts_y_zz[i] * fe_0 + ts_y_yzz[i] * pa_y[i];

        ts_yy_zzz[i] = ts_0_zzz[i] * fe_0 + ts_y_zzz[i] * pa_y[i];
    }

    // Set up 40-50 components of targeted buffer : DF

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

#pragma omp simd aligned(pa_y,          \
                             pa_z,      \
                             ts_y_xxy,  \
                             ts_y_xyy,  \
                             ts_y_yyy,  \
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
                             ts_z_xxx,  \
                             ts_z_xxz,  \
                             ts_z_xyz,  \
                             ts_z_xz,   \
                             ts_z_xzz,  \
                             ts_z_yyz,  \
                             ts_z_yz,   \
                             ts_z_yzz,  \
                             ts_z_zz,   \
                             ts_z_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yz_xxx[i] = ts_z_xxx[i] * pa_y[i];

        ts_yz_xxy[i] = ts_y_xxy[i] * pa_z[i];

        ts_yz_xxz[i] = ts_z_xxz[i] * pa_y[i];

        ts_yz_xyy[i] = ts_y_xyy[i] * pa_z[i];

        ts_yz_xyz[i] = ts_z_xz[i] * fe_0 + ts_z_xyz[i] * pa_y[i];

        ts_yz_xzz[i] = ts_z_xzz[i] * pa_y[i];

        ts_yz_yyy[i] = ts_y_yyy[i] * pa_z[i];

        ts_yz_yyz[i] = 2.0 * ts_z_yz[i] * fe_0 + ts_z_yyz[i] * pa_y[i];

        ts_yz_yzz[i] = ts_z_zz[i] * fe_0 + ts_z_yzz[i] * pa_y[i];

        ts_yz_zzz[i] = ts_z_zzz[i] * pa_y[i];
    }

    // Set up 50-60 components of targeted buffer : DF

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

#pragma omp simd aligned(pa_z,          \
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
                             ts_z_xx,   \
                             ts_z_xxx,  \
                             ts_z_xxy,  \
                             ts_z_xxz,  \
                             ts_z_xy,   \
                             ts_z_xyy,  \
                             ts_z_xyz,  \
                             ts_z_xz,   \
                             ts_z_xzz,  \
                             ts_z_yy,   \
                             ts_z_yyy,  \
                             ts_z_yyz,  \
                             ts_z_yz,   \
                             ts_z_yzz,  \
                             ts_z_zz,   \
                             ts_z_zzz,  \
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

        ts_zz_xxx[i] = ts_0_xxx[i] * fe_0 + ts_z_xxx[i] * pa_z[i];

        ts_zz_xxy[i] = ts_0_xxy[i] * fe_0 + ts_z_xxy[i] * pa_z[i];

        ts_zz_xxz[i] = ts_0_xxz[i] * fe_0 + ts_z_xx[i] * fe_0 + ts_z_xxz[i] * pa_z[i];

        ts_zz_xyy[i] = ts_0_xyy[i] * fe_0 + ts_z_xyy[i] * pa_z[i];

        ts_zz_xyz[i] = ts_0_xyz[i] * fe_0 + ts_z_xy[i] * fe_0 + ts_z_xyz[i] * pa_z[i];

        ts_zz_xzz[i] = ts_0_xzz[i] * fe_0 + 2.0 * ts_z_xz[i] * fe_0 + ts_z_xzz[i] * pa_z[i];

        ts_zz_yyy[i] = ts_0_yyy[i] * fe_0 + ts_z_yyy[i] * pa_z[i];

        ts_zz_yyz[i] = ts_0_yyz[i] * fe_0 + ts_z_yy[i] * fe_0 + ts_z_yyz[i] * pa_z[i];

        ts_zz_yzz[i] = ts_0_yzz[i] * fe_0 + 2.0 * ts_z_yz[i] * fe_0 + ts_z_yzz[i] * pa_z[i];

        ts_zz_zzz[i] = ts_0_zzz[i] * fe_0 + 3.0 * ts_z_zz[i] * fe_0 + ts_z_zzz[i] * pa_z[i];
    }
}

}  // namespace ovlrec
