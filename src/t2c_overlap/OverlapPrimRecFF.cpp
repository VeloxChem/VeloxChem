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

#include "OverlapPrimRecFF.hpp"

namespace ovlrec {  // ovlrec namespace

auto
comp_prim_overlap_ff(CSimdArray<double>&       pbuffer,
                     const size_t              idx_ovl_ff,
                     const size_t              idx_ovl_pf,
                     const size_t              idx_ovl_dd,
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

    // Set up components of auxiliary buffer : DD

    auto ts_xx_xx = pbuffer.data(idx_ovl_dd);

    auto ts_xx_xy = pbuffer.data(idx_ovl_dd + 1);

    auto ts_xx_xz = pbuffer.data(idx_ovl_dd + 2);

    auto ts_xx_yy = pbuffer.data(idx_ovl_dd + 3);

    auto ts_xx_yz = pbuffer.data(idx_ovl_dd + 4);

    auto ts_xx_zz = pbuffer.data(idx_ovl_dd + 5);

    auto ts_yy_xx = pbuffer.data(idx_ovl_dd + 18);

    auto ts_yy_xy = pbuffer.data(idx_ovl_dd + 19);

    auto ts_yy_xz = pbuffer.data(idx_ovl_dd + 20);

    auto ts_yy_yy = pbuffer.data(idx_ovl_dd + 21);

    auto ts_yy_yz = pbuffer.data(idx_ovl_dd + 22);

    auto ts_yy_zz = pbuffer.data(idx_ovl_dd + 23);

    auto ts_yz_yz = pbuffer.data(idx_ovl_dd + 28);

    auto ts_zz_xx = pbuffer.data(idx_ovl_dd + 30);

    auto ts_zz_xy = pbuffer.data(idx_ovl_dd + 31);

    auto ts_zz_xz = pbuffer.data(idx_ovl_dd + 32);

    auto ts_zz_yy = pbuffer.data(idx_ovl_dd + 33);

    auto ts_zz_yz = pbuffer.data(idx_ovl_dd + 34);

    auto ts_zz_zz = pbuffer.data(idx_ovl_dd + 35);

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

    auto ts_xy_xxy = pbuffer.data(idx_ovl_df + 11);

    auto ts_xy_xyy = pbuffer.data(idx_ovl_df + 13);

    auto ts_xy_yyy = pbuffer.data(idx_ovl_df + 16);

    auto ts_xy_yyz = pbuffer.data(idx_ovl_df + 17);

    auto ts_xy_yzz = pbuffer.data(idx_ovl_df + 18);

    auto ts_xz_xxx = pbuffer.data(idx_ovl_df + 20);

    auto ts_xz_xxz = pbuffer.data(idx_ovl_df + 22);

    auto ts_xz_xzz = pbuffer.data(idx_ovl_df + 25);

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

    auto ts_yz_xxz = pbuffer.data(idx_ovl_df + 42);

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

    // Set up 0-10 components of targeted buffer : FF

    auto ts_xxx_xxx = pbuffer.data(idx_ovl_ff);

    auto ts_xxx_xxy = pbuffer.data(idx_ovl_ff + 1);

    auto ts_xxx_xxz = pbuffer.data(idx_ovl_ff + 2);

    auto ts_xxx_xyy = pbuffer.data(idx_ovl_ff + 3);

    auto ts_xxx_xyz = pbuffer.data(idx_ovl_ff + 4);

    auto ts_xxx_xzz = pbuffer.data(idx_ovl_ff + 5);

    auto ts_xxx_yyy = pbuffer.data(idx_ovl_ff + 6);

    auto ts_xxx_yyz = pbuffer.data(idx_ovl_ff + 7);

    auto ts_xxx_yzz = pbuffer.data(idx_ovl_ff + 8);

    auto ts_xxx_zzz = pbuffer.data(idx_ovl_ff + 9);

#pragma omp simd aligned(pa_x,           \
                             ts_x_xxx,   \
                             ts_x_xxy,   \
                             ts_x_xxz,   \
                             ts_x_xyy,   \
                             ts_x_xyz,   \
                             ts_x_xzz,   \
                             ts_x_yyy,   \
                             ts_x_yyz,   \
                             ts_x_yzz,   \
                             ts_x_zzz,   \
                             ts_xx_xx,   \
                             ts_xx_xxx,  \
                             ts_xx_xxy,  \
                             ts_xx_xxz,  \
                             ts_xx_xy,   \
                             ts_xx_xyy,  \
                             ts_xx_xyz,  \
                             ts_xx_xz,   \
                             ts_xx_xzz,  \
                             ts_xx_yy,   \
                             ts_xx_yyy,  \
                             ts_xx_yyz,  \
                             ts_xx_yz,   \
                             ts_xx_yzz,  \
                             ts_xx_zz,   \
                             ts_xx_zzz,  \
                             ts_xxx_xxx, \
                             ts_xxx_xxy, \
                             ts_xxx_xxz, \
                             ts_xxx_xyy, \
                             ts_xxx_xyz, \
                             ts_xxx_xzz, \
                             ts_xxx_yyy, \
                             ts_xxx_yyz, \
                             ts_xxx_yzz, \
                             ts_xxx_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxx_xxx[i] = 2.0 * ts_x_xxx[i] * fe_0 + 3.0 * ts_xx_xx[i] * fe_0 + ts_xx_xxx[i] * pa_x[i];

        ts_xxx_xxy[i] = 2.0 * ts_x_xxy[i] * fe_0 + 2.0 * ts_xx_xy[i] * fe_0 + ts_xx_xxy[i] * pa_x[i];

        ts_xxx_xxz[i] = 2.0 * ts_x_xxz[i] * fe_0 + 2.0 * ts_xx_xz[i] * fe_0 + ts_xx_xxz[i] * pa_x[i];

        ts_xxx_xyy[i] = 2.0 * ts_x_xyy[i] * fe_0 + ts_xx_yy[i] * fe_0 + ts_xx_xyy[i] * pa_x[i];

        ts_xxx_xyz[i] = 2.0 * ts_x_xyz[i] * fe_0 + ts_xx_yz[i] * fe_0 + ts_xx_xyz[i] * pa_x[i];

        ts_xxx_xzz[i] = 2.0 * ts_x_xzz[i] * fe_0 + ts_xx_zz[i] * fe_0 + ts_xx_xzz[i] * pa_x[i];

        ts_xxx_yyy[i] = 2.0 * ts_x_yyy[i] * fe_0 + ts_xx_yyy[i] * pa_x[i];

        ts_xxx_yyz[i] = 2.0 * ts_x_yyz[i] * fe_0 + ts_xx_yyz[i] * pa_x[i];

        ts_xxx_yzz[i] = 2.0 * ts_x_yzz[i] * fe_0 + ts_xx_yzz[i] * pa_x[i];

        ts_xxx_zzz[i] = 2.0 * ts_x_zzz[i] * fe_0 + ts_xx_zzz[i] * pa_x[i];
    }

    // Set up 10-20 components of targeted buffer : FF

    auto ts_xxy_xxx = pbuffer.data(idx_ovl_ff + 10);

    auto ts_xxy_xxy = pbuffer.data(idx_ovl_ff + 11);

    auto ts_xxy_xxz = pbuffer.data(idx_ovl_ff + 12);

    auto ts_xxy_xyy = pbuffer.data(idx_ovl_ff + 13);

    auto ts_xxy_xyz = pbuffer.data(idx_ovl_ff + 14);

    auto ts_xxy_xzz = pbuffer.data(idx_ovl_ff + 15);

    auto ts_xxy_yyy = pbuffer.data(idx_ovl_ff + 16);

    auto ts_xxy_yyz = pbuffer.data(idx_ovl_ff + 17);

    auto ts_xxy_yzz = pbuffer.data(idx_ovl_ff + 18);

    auto ts_xxy_zzz = pbuffer.data(idx_ovl_ff + 19);

#pragma omp simd aligned(pa_x,           \
                             pa_y,       \
                             ts_xx_xx,   \
                             ts_xx_xxx,  \
                             ts_xx_xxy,  \
                             ts_xx_xxz,  \
                             ts_xx_xy,   \
                             ts_xx_xyy,  \
                             ts_xx_xyz,  \
                             ts_xx_xz,   \
                             ts_xx_xzz,  \
                             ts_xx_zzz,  \
                             ts_xxy_xxx, \
                             ts_xxy_xxy, \
                             ts_xxy_xxz, \
                             ts_xxy_xyy, \
                             ts_xxy_xyz, \
                             ts_xxy_xzz, \
                             ts_xxy_yyy, \
                             ts_xxy_yyz, \
                             ts_xxy_yzz, \
                             ts_xxy_zzz, \
                             ts_xy_yyy,  \
                             ts_xy_yyz,  \
                             ts_xy_yzz,  \
                             ts_y_yyy,   \
                             ts_y_yyz,   \
                             ts_y_yzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxy_xxx[i] = ts_xx_xxx[i] * pa_y[i];

        ts_xxy_xxy[i] = ts_xx_xx[i] * fe_0 + ts_xx_xxy[i] * pa_y[i];

        ts_xxy_xxz[i] = ts_xx_xxz[i] * pa_y[i];

        ts_xxy_xyy[i] = 2.0 * ts_xx_xy[i] * fe_0 + ts_xx_xyy[i] * pa_y[i];

        ts_xxy_xyz[i] = ts_xx_xz[i] * fe_0 + ts_xx_xyz[i] * pa_y[i];

        ts_xxy_xzz[i] = ts_xx_xzz[i] * pa_y[i];

        ts_xxy_yyy[i] = ts_y_yyy[i] * fe_0 + ts_xy_yyy[i] * pa_x[i];

        ts_xxy_yyz[i] = ts_y_yyz[i] * fe_0 + ts_xy_yyz[i] * pa_x[i];

        ts_xxy_yzz[i] = ts_y_yzz[i] * fe_0 + ts_xy_yzz[i] * pa_x[i];

        ts_xxy_zzz[i] = ts_xx_zzz[i] * pa_y[i];
    }

    // Set up 20-30 components of targeted buffer : FF

    auto ts_xxz_xxx = pbuffer.data(idx_ovl_ff + 20);

    auto ts_xxz_xxy = pbuffer.data(idx_ovl_ff + 21);

    auto ts_xxz_xxz = pbuffer.data(idx_ovl_ff + 22);

    auto ts_xxz_xyy = pbuffer.data(idx_ovl_ff + 23);

    auto ts_xxz_xyz = pbuffer.data(idx_ovl_ff + 24);

    auto ts_xxz_xzz = pbuffer.data(idx_ovl_ff + 25);

    auto ts_xxz_yyy = pbuffer.data(idx_ovl_ff + 26);

    auto ts_xxz_yyz = pbuffer.data(idx_ovl_ff + 27);

    auto ts_xxz_yzz = pbuffer.data(idx_ovl_ff + 28);

    auto ts_xxz_zzz = pbuffer.data(idx_ovl_ff + 29);

#pragma omp simd aligned(pa_x,           \
                             pa_z,       \
                             ts_xx_xx,   \
                             ts_xx_xxx,  \
                             ts_xx_xxy,  \
                             ts_xx_xxz,  \
                             ts_xx_xy,   \
                             ts_xx_xyy,  \
                             ts_xx_xyz,  \
                             ts_xx_xz,   \
                             ts_xx_xzz,  \
                             ts_xx_yyy,  \
                             ts_xxz_xxx, \
                             ts_xxz_xxy, \
                             ts_xxz_xxz, \
                             ts_xxz_xyy, \
                             ts_xxz_xyz, \
                             ts_xxz_xzz, \
                             ts_xxz_yyy, \
                             ts_xxz_yyz, \
                             ts_xxz_yzz, \
                             ts_xxz_zzz, \
                             ts_xz_yyz,  \
                             ts_xz_yzz,  \
                             ts_xz_zzz,  \
                             ts_z_yyz,   \
                             ts_z_yzz,   \
                             ts_z_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxz_xxx[i] = ts_xx_xxx[i] * pa_z[i];

        ts_xxz_xxy[i] = ts_xx_xxy[i] * pa_z[i];

        ts_xxz_xxz[i] = ts_xx_xx[i] * fe_0 + ts_xx_xxz[i] * pa_z[i];

        ts_xxz_xyy[i] = ts_xx_xyy[i] * pa_z[i];

        ts_xxz_xyz[i] = ts_xx_xy[i] * fe_0 + ts_xx_xyz[i] * pa_z[i];

        ts_xxz_xzz[i] = 2.0 * ts_xx_xz[i] * fe_0 + ts_xx_xzz[i] * pa_z[i];

        ts_xxz_yyy[i] = ts_xx_yyy[i] * pa_z[i];

        ts_xxz_yyz[i] = ts_z_yyz[i] * fe_0 + ts_xz_yyz[i] * pa_x[i];

        ts_xxz_yzz[i] = ts_z_yzz[i] * fe_0 + ts_xz_yzz[i] * pa_x[i];

        ts_xxz_zzz[i] = ts_z_zzz[i] * fe_0 + ts_xz_zzz[i] * pa_x[i];
    }

    // Set up 30-40 components of targeted buffer : FF

    auto ts_xyy_xxx = pbuffer.data(idx_ovl_ff + 30);

    auto ts_xyy_xxy = pbuffer.data(idx_ovl_ff + 31);

    auto ts_xyy_xxz = pbuffer.data(idx_ovl_ff + 32);

    auto ts_xyy_xyy = pbuffer.data(idx_ovl_ff + 33);

    auto ts_xyy_xyz = pbuffer.data(idx_ovl_ff + 34);

    auto ts_xyy_xzz = pbuffer.data(idx_ovl_ff + 35);

    auto ts_xyy_yyy = pbuffer.data(idx_ovl_ff + 36);

    auto ts_xyy_yyz = pbuffer.data(idx_ovl_ff + 37);

    auto ts_xyy_yzz = pbuffer.data(idx_ovl_ff + 38);

    auto ts_xyy_zzz = pbuffer.data(idx_ovl_ff + 39);

#pragma omp simd aligned(pa_x,           \
                             ts_xyy_xxx, \
                             ts_xyy_xxy, \
                             ts_xyy_xxz, \
                             ts_xyy_xyy, \
                             ts_xyy_xyz, \
                             ts_xyy_xzz, \
                             ts_xyy_yyy, \
                             ts_xyy_yyz, \
                             ts_xyy_yzz, \
                             ts_xyy_zzz, \
                             ts_yy_xx,   \
                             ts_yy_xxx,  \
                             ts_yy_xxy,  \
                             ts_yy_xxz,  \
                             ts_yy_xy,   \
                             ts_yy_xyy,  \
                             ts_yy_xyz,  \
                             ts_yy_xz,   \
                             ts_yy_xzz,  \
                             ts_yy_yy,   \
                             ts_yy_yyy,  \
                             ts_yy_yyz,  \
                             ts_yy_yz,   \
                             ts_yy_yzz,  \
                             ts_yy_zz,   \
                             ts_yy_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyy_xxx[i] = 3.0 * ts_yy_xx[i] * fe_0 + ts_yy_xxx[i] * pa_x[i];

        ts_xyy_xxy[i] = 2.0 * ts_yy_xy[i] * fe_0 + ts_yy_xxy[i] * pa_x[i];

        ts_xyy_xxz[i] = 2.0 * ts_yy_xz[i] * fe_0 + ts_yy_xxz[i] * pa_x[i];

        ts_xyy_xyy[i] = ts_yy_yy[i] * fe_0 + ts_yy_xyy[i] * pa_x[i];

        ts_xyy_xyz[i] = ts_yy_yz[i] * fe_0 + ts_yy_xyz[i] * pa_x[i];

        ts_xyy_xzz[i] = ts_yy_zz[i] * fe_0 + ts_yy_xzz[i] * pa_x[i];

        ts_xyy_yyy[i] = ts_yy_yyy[i] * pa_x[i];

        ts_xyy_yyz[i] = ts_yy_yyz[i] * pa_x[i];

        ts_xyy_yzz[i] = ts_yy_yzz[i] * pa_x[i];

        ts_xyy_zzz[i] = ts_yy_zzz[i] * pa_x[i];
    }

    // Set up 40-50 components of targeted buffer : FF

    auto ts_xyz_xxx = pbuffer.data(idx_ovl_ff + 40);

    auto ts_xyz_xxy = pbuffer.data(idx_ovl_ff + 41);

    auto ts_xyz_xxz = pbuffer.data(idx_ovl_ff + 42);

    auto ts_xyz_xyy = pbuffer.data(idx_ovl_ff + 43);

    auto ts_xyz_xyz = pbuffer.data(idx_ovl_ff + 44);

    auto ts_xyz_xzz = pbuffer.data(idx_ovl_ff + 45);

    auto ts_xyz_yyy = pbuffer.data(idx_ovl_ff + 46);

    auto ts_xyz_yyz = pbuffer.data(idx_ovl_ff + 47);

    auto ts_xyz_yzz = pbuffer.data(idx_ovl_ff + 48);

    auto ts_xyz_zzz = pbuffer.data(idx_ovl_ff + 49);

#pragma omp simd aligned(pa_x,           \
                             pa_y,       \
                             pa_z,       \
                             ts_xy_xxy,  \
                             ts_xy_xyy,  \
                             ts_xyz_xxx, \
                             ts_xyz_xxy, \
                             ts_xyz_xxz, \
                             ts_xyz_xyy, \
                             ts_xyz_xyz, \
                             ts_xyz_xzz, \
                             ts_xyz_yyy, \
                             ts_xyz_yyz, \
                             ts_xyz_yzz, \
                             ts_xyz_zzz, \
                             ts_xz_xxx,  \
                             ts_xz_xxz,  \
                             ts_xz_xzz,  \
                             ts_yz_xyz,  \
                             ts_yz_yyy,  \
                             ts_yz_yyz,  \
                             ts_yz_yz,   \
                             ts_yz_yzz,  \
                             ts_yz_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyz_xxx[i] = ts_xz_xxx[i] * pa_y[i];

        ts_xyz_xxy[i] = ts_xy_xxy[i] * pa_z[i];

        ts_xyz_xxz[i] = ts_xz_xxz[i] * pa_y[i];

        ts_xyz_xyy[i] = ts_xy_xyy[i] * pa_z[i];

        ts_xyz_xyz[i] = ts_yz_yz[i] * fe_0 + ts_yz_xyz[i] * pa_x[i];

        ts_xyz_xzz[i] = ts_xz_xzz[i] * pa_y[i];

        ts_xyz_yyy[i] = ts_yz_yyy[i] * pa_x[i];

        ts_xyz_yyz[i] = ts_yz_yyz[i] * pa_x[i];

        ts_xyz_yzz[i] = ts_yz_yzz[i] * pa_x[i];

        ts_xyz_zzz[i] = ts_yz_zzz[i] * pa_x[i];
    }

    // Set up 50-60 components of targeted buffer : FF

    auto ts_xzz_xxx = pbuffer.data(idx_ovl_ff + 50);

    auto ts_xzz_xxy = pbuffer.data(idx_ovl_ff + 51);

    auto ts_xzz_xxz = pbuffer.data(idx_ovl_ff + 52);

    auto ts_xzz_xyy = pbuffer.data(idx_ovl_ff + 53);

    auto ts_xzz_xyz = pbuffer.data(idx_ovl_ff + 54);

    auto ts_xzz_xzz = pbuffer.data(idx_ovl_ff + 55);

    auto ts_xzz_yyy = pbuffer.data(idx_ovl_ff + 56);

    auto ts_xzz_yyz = pbuffer.data(idx_ovl_ff + 57);

    auto ts_xzz_yzz = pbuffer.data(idx_ovl_ff + 58);

    auto ts_xzz_zzz = pbuffer.data(idx_ovl_ff + 59);

#pragma omp simd aligned(pa_x,           \
                             ts_xzz_xxx, \
                             ts_xzz_xxy, \
                             ts_xzz_xxz, \
                             ts_xzz_xyy, \
                             ts_xzz_xyz, \
                             ts_xzz_xzz, \
                             ts_xzz_yyy, \
                             ts_xzz_yyz, \
                             ts_xzz_yzz, \
                             ts_xzz_zzz, \
                             ts_zz_xx,   \
                             ts_zz_xxx,  \
                             ts_zz_xxy,  \
                             ts_zz_xxz,  \
                             ts_zz_xy,   \
                             ts_zz_xyy,  \
                             ts_zz_xyz,  \
                             ts_zz_xz,   \
                             ts_zz_xzz,  \
                             ts_zz_yy,   \
                             ts_zz_yyy,  \
                             ts_zz_yyz,  \
                             ts_zz_yz,   \
                             ts_zz_yzz,  \
                             ts_zz_zz,   \
                             ts_zz_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xzz_xxx[i] = 3.0 * ts_zz_xx[i] * fe_0 + ts_zz_xxx[i] * pa_x[i];

        ts_xzz_xxy[i] = 2.0 * ts_zz_xy[i] * fe_0 + ts_zz_xxy[i] * pa_x[i];

        ts_xzz_xxz[i] = 2.0 * ts_zz_xz[i] * fe_0 + ts_zz_xxz[i] * pa_x[i];

        ts_xzz_xyy[i] = ts_zz_yy[i] * fe_0 + ts_zz_xyy[i] * pa_x[i];

        ts_xzz_xyz[i] = ts_zz_yz[i] * fe_0 + ts_zz_xyz[i] * pa_x[i];

        ts_xzz_xzz[i] = ts_zz_zz[i] * fe_0 + ts_zz_xzz[i] * pa_x[i];

        ts_xzz_yyy[i] = ts_zz_yyy[i] * pa_x[i];

        ts_xzz_yyz[i] = ts_zz_yyz[i] * pa_x[i];

        ts_xzz_yzz[i] = ts_zz_yzz[i] * pa_x[i];

        ts_xzz_zzz[i] = ts_zz_zzz[i] * pa_x[i];
    }

    // Set up 60-70 components of targeted buffer : FF

    auto ts_yyy_xxx = pbuffer.data(idx_ovl_ff + 60);

    auto ts_yyy_xxy = pbuffer.data(idx_ovl_ff + 61);

    auto ts_yyy_xxz = pbuffer.data(idx_ovl_ff + 62);

    auto ts_yyy_xyy = pbuffer.data(idx_ovl_ff + 63);

    auto ts_yyy_xyz = pbuffer.data(idx_ovl_ff + 64);

    auto ts_yyy_xzz = pbuffer.data(idx_ovl_ff + 65);

    auto ts_yyy_yyy = pbuffer.data(idx_ovl_ff + 66);

    auto ts_yyy_yyz = pbuffer.data(idx_ovl_ff + 67);

    auto ts_yyy_yzz = pbuffer.data(idx_ovl_ff + 68);

    auto ts_yyy_zzz = pbuffer.data(idx_ovl_ff + 69);

#pragma omp simd aligned(pa_y,           \
                             ts_y_xxx,   \
                             ts_y_xxy,   \
                             ts_y_xxz,   \
                             ts_y_xyy,   \
                             ts_y_xyz,   \
                             ts_y_xzz,   \
                             ts_y_yyy,   \
                             ts_y_yyz,   \
                             ts_y_yzz,   \
                             ts_y_zzz,   \
                             ts_yy_xx,   \
                             ts_yy_xxx,  \
                             ts_yy_xxy,  \
                             ts_yy_xxz,  \
                             ts_yy_xy,   \
                             ts_yy_xyy,  \
                             ts_yy_xyz,  \
                             ts_yy_xz,   \
                             ts_yy_xzz,  \
                             ts_yy_yy,   \
                             ts_yy_yyy,  \
                             ts_yy_yyz,  \
                             ts_yy_yz,   \
                             ts_yy_yzz,  \
                             ts_yy_zz,   \
                             ts_yy_zzz,  \
                             ts_yyy_xxx, \
                             ts_yyy_xxy, \
                             ts_yyy_xxz, \
                             ts_yyy_xyy, \
                             ts_yyy_xyz, \
                             ts_yyy_xzz, \
                             ts_yyy_yyy, \
                             ts_yyy_yyz, \
                             ts_yyy_yzz, \
                             ts_yyy_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyy_xxx[i] = 2.0 * ts_y_xxx[i] * fe_0 + ts_yy_xxx[i] * pa_y[i];

        ts_yyy_xxy[i] = 2.0 * ts_y_xxy[i] * fe_0 + ts_yy_xx[i] * fe_0 + ts_yy_xxy[i] * pa_y[i];

        ts_yyy_xxz[i] = 2.0 * ts_y_xxz[i] * fe_0 + ts_yy_xxz[i] * pa_y[i];

        ts_yyy_xyy[i] = 2.0 * ts_y_xyy[i] * fe_0 + 2.0 * ts_yy_xy[i] * fe_0 + ts_yy_xyy[i] * pa_y[i];

        ts_yyy_xyz[i] = 2.0 * ts_y_xyz[i] * fe_0 + ts_yy_xz[i] * fe_0 + ts_yy_xyz[i] * pa_y[i];

        ts_yyy_xzz[i] = 2.0 * ts_y_xzz[i] * fe_0 + ts_yy_xzz[i] * pa_y[i];

        ts_yyy_yyy[i] = 2.0 * ts_y_yyy[i] * fe_0 + 3.0 * ts_yy_yy[i] * fe_0 + ts_yy_yyy[i] * pa_y[i];

        ts_yyy_yyz[i] = 2.0 * ts_y_yyz[i] * fe_0 + 2.0 * ts_yy_yz[i] * fe_0 + ts_yy_yyz[i] * pa_y[i];

        ts_yyy_yzz[i] = 2.0 * ts_y_yzz[i] * fe_0 + ts_yy_zz[i] * fe_0 + ts_yy_yzz[i] * pa_y[i];

        ts_yyy_zzz[i] = 2.0 * ts_y_zzz[i] * fe_0 + ts_yy_zzz[i] * pa_y[i];
    }

    // Set up 70-80 components of targeted buffer : FF

    auto ts_yyz_xxx = pbuffer.data(idx_ovl_ff + 70);

    auto ts_yyz_xxy = pbuffer.data(idx_ovl_ff + 71);

    auto ts_yyz_xxz = pbuffer.data(idx_ovl_ff + 72);

    auto ts_yyz_xyy = pbuffer.data(idx_ovl_ff + 73);

    auto ts_yyz_xyz = pbuffer.data(idx_ovl_ff + 74);

    auto ts_yyz_xzz = pbuffer.data(idx_ovl_ff + 75);

    auto ts_yyz_yyy = pbuffer.data(idx_ovl_ff + 76);

    auto ts_yyz_yyz = pbuffer.data(idx_ovl_ff + 77);

    auto ts_yyz_yzz = pbuffer.data(idx_ovl_ff + 78);

    auto ts_yyz_zzz = pbuffer.data(idx_ovl_ff + 79);

#pragma omp simd aligned(pa_y,           \
                             pa_z,       \
                             ts_yy_xxx,  \
                             ts_yy_xxy,  \
                             ts_yy_xy,   \
                             ts_yy_xyy,  \
                             ts_yy_xyz,  \
                             ts_yy_yy,   \
                             ts_yy_yyy,  \
                             ts_yy_yyz,  \
                             ts_yy_yz,   \
                             ts_yy_yzz,  \
                             ts_yyz_xxx, \
                             ts_yyz_xxy, \
                             ts_yyz_xxz, \
                             ts_yyz_xyy, \
                             ts_yyz_xyz, \
                             ts_yyz_xzz, \
                             ts_yyz_yyy, \
                             ts_yyz_yyz, \
                             ts_yyz_yzz, \
                             ts_yyz_zzz, \
                             ts_yz_xxz,  \
                             ts_yz_xzz,  \
                             ts_yz_zzz,  \
                             ts_z_xxz,   \
                             ts_z_xzz,   \
                             ts_z_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyz_xxx[i] = ts_yy_xxx[i] * pa_z[i];

        ts_yyz_xxy[i] = ts_yy_xxy[i] * pa_z[i];

        ts_yyz_xxz[i] = ts_z_xxz[i] * fe_0 + ts_yz_xxz[i] * pa_y[i];

        ts_yyz_xyy[i] = ts_yy_xyy[i] * pa_z[i];

        ts_yyz_xyz[i] = ts_yy_xy[i] * fe_0 + ts_yy_xyz[i] * pa_z[i];

        ts_yyz_xzz[i] = ts_z_xzz[i] * fe_0 + ts_yz_xzz[i] * pa_y[i];

        ts_yyz_yyy[i] = ts_yy_yyy[i] * pa_z[i];

        ts_yyz_yyz[i] = ts_yy_yy[i] * fe_0 + ts_yy_yyz[i] * pa_z[i];

        ts_yyz_yzz[i] = 2.0 * ts_yy_yz[i] * fe_0 + ts_yy_yzz[i] * pa_z[i];

        ts_yyz_zzz[i] = ts_z_zzz[i] * fe_0 + ts_yz_zzz[i] * pa_y[i];
    }

    // Set up 80-90 components of targeted buffer : FF

    auto ts_yzz_xxx = pbuffer.data(idx_ovl_ff + 80);

    auto ts_yzz_xxy = pbuffer.data(idx_ovl_ff + 81);

    auto ts_yzz_xxz = pbuffer.data(idx_ovl_ff + 82);

    auto ts_yzz_xyy = pbuffer.data(idx_ovl_ff + 83);

    auto ts_yzz_xyz = pbuffer.data(idx_ovl_ff + 84);

    auto ts_yzz_xzz = pbuffer.data(idx_ovl_ff + 85);

    auto ts_yzz_yyy = pbuffer.data(idx_ovl_ff + 86);

    auto ts_yzz_yyz = pbuffer.data(idx_ovl_ff + 87);

    auto ts_yzz_yzz = pbuffer.data(idx_ovl_ff + 88);

    auto ts_yzz_zzz = pbuffer.data(idx_ovl_ff + 89);

#pragma omp simd aligned(pa_y,           \
                             ts_yzz_xxx, \
                             ts_yzz_xxy, \
                             ts_yzz_xxz, \
                             ts_yzz_xyy, \
                             ts_yzz_xyz, \
                             ts_yzz_xzz, \
                             ts_yzz_yyy, \
                             ts_yzz_yyz, \
                             ts_yzz_yzz, \
                             ts_yzz_zzz, \
                             ts_zz_xx,   \
                             ts_zz_xxx,  \
                             ts_zz_xxy,  \
                             ts_zz_xxz,  \
                             ts_zz_xy,   \
                             ts_zz_xyy,  \
                             ts_zz_xyz,  \
                             ts_zz_xz,   \
                             ts_zz_xzz,  \
                             ts_zz_yy,   \
                             ts_zz_yyy,  \
                             ts_zz_yyz,  \
                             ts_zz_yz,   \
                             ts_zz_yzz,  \
                             ts_zz_zz,   \
                             ts_zz_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yzz_xxx[i] = ts_zz_xxx[i] * pa_y[i];

        ts_yzz_xxy[i] = ts_zz_xx[i] * fe_0 + ts_zz_xxy[i] * pa_y[i];

        ts_yzz_xxz[i] = ts_zz_xxz[i] * pa_y[i];

        ts_yzz_xyy[i] = 2.0 * ts_zz_xy[i] * fe_0 + ts_zz_xyy[i] * pa_y[i];

        ts_yzz_xyz[i] = ts_zz_xz[i] * fe_0 + ts_zz_xyz[i] * pa_y[i];

        ts_yzz_xzz[i] = ts_zz_xzz[i] * pa_y[i];

        ts_yzz_yyy[i] = 3.0 * ts_zz_yy[i] * fe_0 + ts_zz_yyy[i] * pa_y[i];

        ts_yzz_yyz[i] = 2.0 * ts_zz_yz[i] * fe_0 + ts_zz_yyz[i] * pa_y[i];

        ts_yzz_yzz[i] = ts_zz_zz[i] * fe_0 + ts_zz_yzz[i] * pa_y[i];

        ts_yzz_zzz[i] = ts_zz_zzz[i] * pa_y[i];
    }

    // Set up 90-100 components of targeted buffer : FF

    auto ts_zzz_xxx = pbuffer.data(idx_ovl_ff + 90);

    auto ts_zzz_xxy = pbuffer.data(idx_ovl_ff + 91);

    auto ts_zzz_xxz = pbuffer.data(idx_ovl_ff + 92);

    auto ts_zzz_xyy = pbuffer.data(idx_ovl_ff + 93);

    auto ts_zzz_xyz = pbuffer.data(idx_ovl_ff + 94);

    auto ts_zzz_xzz = pbuffer.data(idx_ovl_ff + 95);

    auto ts_zzz_yyy = pbuffer.data(idx_ovl_ff + 96);

    auto ts_zzz_yyz = pbuffer.data(idx_ovl_ff + 97);

    auto ts_zzz_yzz = pbuffer.data(idx_ovl_ff + 98);

    auto ts_zzz_zzz = pbuffer.data(idx_ovl_ff + 99);

#pragma omp simd aligned(pa_z,           \
                             ts_z_xxx,   \
                             ts_z_xxy,   \
                             ts_z_xxz,   \
                             ts_z_xyy,   \
                             ts_z_xyz,   \
                             ts_z_xzz,   \
                             ts_z_yyy,   \
                             ts_z_yyz,   \
                             ts_z_yzz,   \
                             ts_z_zzz,   \
                             ts_zz_xx,   \
                             ts_zz_xxx,  \
                             ts_zz_xxy,  \
                             ts_zz_xxz,  \
                             ts_zz_xy,   \
                             ts_zz_xyy,  \
                             ts_zz_xyz,  \
                             ts_zz_xz,   \
                             ts_zz_xzz,  \
                             ts_zz_yy,   \
                             ts_zz_yyy,  \
                             ts_zz_yyz,  \
                             ts_zz_yz,   \
                             ts_zz_yzz,  \
                             ts_zz_zz,   \
                             ts_zz_zzz,  \
                             ts_zzz_xxx, \
                             ts_zzz_xxy, \
                             ts_zzz_xxz, \
                             ts_zzz_xyy, \
                             ts_zzz_xyz, \
                             ts_zzz_xzz, \
                             ts_zzz_yyy, \
                             ts_zzz_yyz, \
                             ts_zzz_yzz, \
                             ts_zzz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_zzz_xxx[i] = 2.0 * ts_z_xxx[i] * fe_0 + ts_zz_xxx[i] * pa_z[i];

        ts_zzz_xxy[i] = 2.0 * ts_z_xxy[i] * fe_0 + ts_zz_xxy[i] * pa_z[i];

        ts_zzz_xxz[i] = 2.0 * ts_z_xxz[i] * fe_0 + ts_zz_xx[i] * fe_0 + ts_zz_xxz[i] * pa_z[i];

        ts_zzz_xyy[i] = 2.0 * ts_z_xyy[i] * fe_0 + ts_zz_xyy[i] * pa_z[i];

        ts_zzz_xyz[i] = 2.0 * ts_z_xyz[i] * fe_0 + ts_zz_xy[i] * fe_0 + ts_zz_xyz[i] * pa_z[i];

        ts_zzz_xzz[i] = 2.0 * ts_z_xzz[i] * fe_0 + 2.0 * ts_zz_xz[i] * fe_0 + ts_zz_xzz[i] * pa_z[i];

        ts_zzz_yyy[i] = 2.0 * ts_z_yyy[i] * fe_0 + ts_zz_yyy[i] * pa_z[i];

        ts_zzz_yyz[i] = 2.0 * ts_z_yyz[i] * fe_0 + ts_zz_yy[i] * fe_0 + ts_zz_yyz[i] * pa_z[i];

        ts_zzz_yzz[i] = 2.0 * ts_z_yzz[i] * fe_0 + 2.0 * ts_zz_yz[i] * fe_0 + ts_zz_yzz[i] * pa_z[i];

        ts_zzz_zzz[i] = 2.0 * ts_z_zzz[i] * fe_0 + 3.0 * ts_zz_zz[i] * fe_0 + ts_zz_zzz[i] * pa_z[i];
    }
}

}  // namespace ovlrec
