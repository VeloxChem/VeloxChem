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

#include "OverlapPrimRecFD.hpp"

namespace ovlrec {  // ovlrec namespace

auto
comp_prim_overlap_fd(CSimdArray<double>&       pbuffer,
                     const size_t              idx_ovl_fd,
                     const size_t              idx_ovl_pd,
                     const size_t              idx_ovl_dp,
                     const size_t              idx_ovl_dd,
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

    // Set up components of auxiliary buffer : DP

    auto ts_xx_x = pbuffer.data(idx_ovl_dp);

    auto ts_xx_y = pbuffer.data(idx_ovl_dp + 1);

    auto ts_xx_z = pbuffer.data(idx_ovl_dp + 2);

    auto ts_yy_x = pbuffer.data(idx_ovl_dp + 9);

    auto ts_yy_y = pbuffer.data(idx_ovl_dp + 10);

    auto ts_yy_z = pbuffer.data(idx_ovl_dp + 11);

    auto ts_zz_x = pbuffer.data(idx_ovl_dp + 15);

    auto ts_zz_y = pbuffer.data(idx_ovl_dp + 16);

    auto ts_zz_z = pbuffer.data(idx_ovl_dp + 17);

    // Set up components of auxiliary buffer : DD

    auto ts_xx_xx = pbuffer.data(idx_ovl_dd);

    auto ts_xx_xy = pbuffer.data(idx_ovl_dd + 1);

    auto ts_xx_xz = pbuffer.data(idx_ovl_dd + 2);

    auto ts_xx_yy = pbuffer.data(idx_ovl_dd + 3);

    auto ts_xx_yz = pbuffer.data(idx_ovl_dd + 4);

    auto ts_xx_zz = pbuffer.data(idx_ovl_dd + 5);

    auto ts_xy_xy = pbuffer.data(idx_ovl_dd + 7);

    auto ts_xy_yy = pbuffer.data(idx_ovl_dd + 9);

    auto ts_xy_yz = pbuffer.data(idx_ovl_dd + 10);

    auto ts_xz_xx = pbuffer.data(idx_ovl_dd + 12);

    auto ts_xz_xz = pbuffer.data(idx_ovl_dd + 14);

    auto ts_xz_yz = pbuffer.data(idx_ovl_dd + 16);

    auto ts_xz_zz = pbuffer.data(idx_ovl_dd + 17);

    auto ts_yy_xx = pbuffer.data(idx_ovl_dd + 18);

    auto ts_yy_xy = pbuffer.data(idx_ovl_dd + 19);

    auto ts_yy_xz = pbuffer.data(idx_ovl_dd + 20);

    auto ts_yy_yy = pbuffer.data(idx_ovl_dd + 21);

    auto ts_yy_yz = pbuffer.data(idx_ovl_dd + 22);

    auto ts_yy_zz = pbuffer.data(idx_ovl_dd + 23);

    auto ts_yz_xz = pbuffer.data(idx_ovl_dd + 26);

    auto ts_yz_yy = pbuffer.data(idx_ovl_dd + 27);

    auto ts_yz_yz = pbuffer.data(idx_ovl_dd + 28);

    auto ts_yz_zz = pbuffer.data(idx_ovl_dd + 29);

    auto ts_zz_xx = pbuffer.data(idx_ovl_dd + 30);

    auto ts_zz_xy = pbuffer.data(idx_ovl_dd + 31);

    auto ts_zz_xz = pbuffer.data(idx_ovl_dd + 32);

    auto ts_zz_yy = pbuffer.data(idx_ovl_dd + 33);

    auto ts_zz_yz = pbuffer.data(idx_ovl_dd + 34);

    auto ts_zz_zz = pbuffer.data(idx_ovl_dd + 35);

    // Set up 0-6 components of targeted buffer : FD

    auto ts_xxx_xx = pbuffer.data(idx_ovl_fd);

    auto ts_xxx_xy = pbuffer.data(idx_ovl_fd + 1);

    auto ts_xxx_xz = pbuffer.data(idx_ovl_fd + 2);

    auto ts_xxx_yy = pbuffer.data(idx_ovl_fd + 3);

    auto ts_xxx_yz = pbuffer.data(idx_ovl_fd + 4);

    auto ts_xxx_zz = pbuffer.data(idx_ovl_fd + 5);

#pragma omp simd aligned(pa_x,          \
                             ts_x_xx,   \
                             ts_x_xy,   \
                             ts_x_xz,   \
                             ts_x_yy,   \
                             ts_x_yz,   \
                             ts_x_zz,   \
                             ts_xx_x,   \
                             ts_xx_xx,  \
                             ts_xx_xy,  \
                             ts_xx_xz,  \
                             ts_xx_y,   \
                             ts_xx_yy,  \
                             ts_xx_yz,  \
                             ts_xx_z,   \
                             ts_xx_zz,  \
                             ts_xxx_xx, \
                             ts_xxx_xy, \
                             ts_xxx_xz, \
                             ts_xxx_yy, \
                             ts_xxx_yz, \
                             ts_xxx_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxx_xx[i] = 2.0 * ts_x_xx[i] * fe_0 + 2.0 * ts_xx_x[i] * fe_0 + ts_xx_xx[i] * pa_x[i];

        ts_xxx_xy[i] = 2.0 * ts_x_xy[i] * fe_0 + ts_xx_y[i] * fe_0 + ts_xx_xy[i] * pa_x[i];

        ts_xxx_xz[i] = 2.0 * ts_x_xz[i] * fe_0 + ts_xx_z[i] * fe_0 + ts_xx_xz[i] * pa_x[i];

        ts_xxx_yy[i] = 2.0 * ts_x_yy[i] * fe_0 + ts_xx_yy[i] * pa_x[i];

        ts_xxx_yz[i] = 2.0 * ts_x_yz[i] * fe_0 + ts_xx_yz[i] * pa_x[i];

        ts_xxx_zz[i] = 2.0 * ts_x_zz[i] * fe_0 + ts_xx_zz[i] * pa_x[i];
    }

    // Set up 6-12 components of targeted buffer : FD

    auto ts_xxy_xx = pbuffer.data(idx_ovl_fd + 6);

    auto ts_xxy_xy = pbuffer.data(idx_ovl_fd + 7);

    auto ts_xxy_xz = pbuffer.data(idx_ovl_fd + 8);

    auto ts_xxy_yy = pbuffer.data(idx_ovl_fd + 9);

    auto ts_xxy_yz = pbuffer.data(idx_ovl_fd + 10);

    auto ts_xxy_zz = pbuffer.data(idx_ovl_fd + 11);

#pragma omp simd aligned(pa_x,          \
                             pa_y,      \
                             ts_xx_x,   \
                             ts_xx_xx,  \
                             ts_xx_xy,  \
                             ts_xx_xz,  \
                             ts_xx_zz,  \
                             ts_xxy_xx, \
                             ts_xxy_xy, \
                             ts_xxy_xz, \
                             ts_xxy_yy, \
                             ts_xxy_yz, \
                             ts_xxy_zz, \
                             ts_xy_yy,  \
                             ts_xy_yz,  \
                             ts_y_yy,   \
                             ts_y_yz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxy_xx[i] = ts_xx_xx[i] * pa_y[i];

        ts_xxy_xy[i] = ts_xx_x[i] * fe_0 + ts_xx_xy[i] * pa_y[i];

        ts_xxy_xz[i] = ts_xx_xz[i] * pa_y[i];

        ts_xxy_yy[i] = ts_y_yy[i] * fe_0 + ts_xy_yy[i] * pa_x[i];

        ts_xxy_yz[i] = ts_y_yz[i] * fe_0 + ts_xy_yz[i] * pa_x[i];

        ts_xxy_zz[i] = ts_xx_zz[i] * pa_y[i];
    }

    // Set up 12-18 components of targeted buffer : FD

    auto ts_xxz_xx = pbuffer.data(idx_ovl_fd + 12);

    auto ts_xxz_xy = pbuffer.data(idx_ovl_fd + 13);

    auto ts_xxz_xz = pbuffer.data(idx_ovl_fd + 14);

    auto ts_xxz_yy = pbuffer.data(idx_ovl_fd + 15);

    auto ts_xxz_yz = pbuffer.data(idx_ovl_fd + 16);

    auto ts_xxz_zz = pbuffer.data(idx_ovl_fd + 17);

#pragma omp simd aligned(pa_x,          \
                             pa_z,      \
                             ts_xx_x,   \
                             ts_xx_xx,  \
                             ts_xx_xy,  \
                             ts_xx_xz,  \
                             ts_xx_yy,  \
                             ts_xxz_xx, \
                             ts_xxz_xy, \
                             ts_xxz_xz, \
                             ts_xxz_yy, \
                             ts_xxz_yz, \
                             ts_xxz_zz, \
                             ts_xz_yz,  \
                             ts_xz_zz,  \
                             ts_z_yz,   \
                             ts_z_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxz_xx[i] = ts_xx_xx[i] * pa_z[i];

        ts_xxz_xy[i] = ts_xx_xy[i] * pa_z[i];

        ts_xxz_xz[i] = ts_xx_x[i] * fe_0 + ts_xx_xz[i] * pa_z[i];

        ts_xxz_yy[i] = ts_xx_yy[i] * pa_z[i];

        ts_xxz_yz[i] = ts_z_yz[i] * fe_0 + ts_xz_yz[i] * pa_x[i];

        ts_xxz_zz[i] = ts_z_zz[i] * fe_0 + ts_xz_zz[i] * pa_x[i];
    }

    // Set up 18-24 components of targeted buffer : FD

    auto ts_xyy_xx = pbuffer.data(idx_ovl_fd + 18);

    auto ts_xyy_xy = pbuffer.data(idx_ovl_fd + 19);

    auto ts_xyy_xz = pbuffer.data(idx_ovl_fd + 20);

    auto ts_xyy_yy = pbuffer.data(idx_ovl_fd + 21);

    auto ts_xyy_yz = pbuffer.data(idx_ovl_fd + 22);

    auto ts_xyy_zz = pbuffer.data(idx_ovl_fd + 23);

#pragma omp simd aligned(pa_x,          \
                             ts_xyy_xx, \
                             ts_xyy_xy, \
                             ts_xyy_xz, \
                             ts_xyy_yy, \
                             ts_xyy_yz, \
                             ts_xyy_zz, \
                             ts_yy_x,   \
                             ts_yy_xx,  \
                             ts_yy_xy,  \
                             ts_yy_xz,  \
                             ts_yy_y,   \
                             ts_yy_yy,  \
                             ts_yy_yz,  \
                             ts_yy_z,   \
                             ts_yy_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyy_xx[i] = 2.0 * ts_yy_x[i] * fe_0 + ts_yy_xx[i] * pa_x[i];

        ts_xyy_xy[i] = ts_yy_y[i] * fe_0 + ts_yy_xy[i] * pa_x[i];

        ts_xyy_xz[i] = ts_yy_z[i] * fe_0 + ts_yy_xz[i] * pa_x[i];

        ts_xyy_yy[i] = ts_yy_yy[i] * pa_x[i];

        ts_xyy_yz[i] = ts_yy_yz[i] * pa_x[i];

        ts_xyy_zz[i] = ts_yy_zz[i] * pa_x[i];
    }

    // Set up 24-30 components of targeted buffer : FD

    auto ts_xyz_xx = pbuffer.data(idx_ovl_fd + 24);

    auto ts_xyz_xy = pbuffer.data(idx_ovl_fd + 25);

    auto ts_xyz_xz = pbuffer.data(idx_ovl_fd + 26);

    auto ts_xyz_yy = pbuffer.data(idx_ovl_fd + 27);

    auto ts_xyz_yz = pbuffer.data(idx_ovl_fd + 28);

    auto ts_xyz_zz = pbuffer.data(idx_ovl_fd + 29);

#pragma omp simd aligned(pa_x,          \
                             pa_y,      \
                             pa_z,      \
                             ts_xy_xy,  \
                             ts_xyz_xx, \
                             ts_xyz_xy, \
                             ts_xyz_xz, \
                             ts_xyz_yy, \
                             ts_xyz_yz, \
                             ts_xyz_zz, \
                             ts_xz_xx,  \
                             ts_xz_xz,  \
                             ts_yz_yy,  \
                             ts_yz_yz,  \
                             ts_yz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ts_xyz_xx[i] = ts_xz_xx[i] * pa_y[i];

        ts_xyz_xy[i] = ts_xy_xy[i] * pa_z[i];

        ts_xyz_xz[i] = ts_xz_xz[i] * pa_y[i];

        ts_xyz_yy[i] = ts_yz_yy[i] * pa_x[i];

        ts_xyz_yz[i] = ts_yz_yz[i] * pa_x[i];

        ts_xyz_zz[i] = ts_yz_zz[i] * pa_x[i];
    }

    // Set up 30-36 components of targeted buffer : FD

    auto ts_xzz_xx = pbuffer.data(idx_ovl_fd + 30);

    auto ts_xzz_xy = pbuffer.data(idx_ovl_fd + 31);

    auto ts_xzz_xz = pbuffer.data(idx_ovl_fd + 32);

    auto ts_xzz_yy = pbuffer.data(idx_ovl_fd + 33);

    auto ts_xzz_yz = pbuffer.data(idx_ovl_fd + 34);

    auto ts_xzz_zz = pbuffer.data(idx_ovl_fd + 35);

#pragma omp simd aligned(pa_x,          \
                             ts_xzz_xx, \
                             ts_xzz_xy, \
                             ts_xzz_xz, \
                             ts_xzz_yy, \
                             ts_xzz_yz, \
                             ts_xzz_zz, \
                             ts_zz_x,   \
                             ts_zz_xx,  \
                             ts_zz_xy,  \
                             ts_zz_xz,  \
                             ts_zz_y,   \
                             ts_zz_yy,  \
                             ts_zz_yz,  \
                             ts_zz_z,   \
                             ts_zz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xzz_xx[i] = 2.0 * ts_zz_x[i] * fe_0 + ts_zz_xx[i] * pa_x[i];

        ts_xzz_xy[i] = ts_zz_y[i] * fe_0 + ts_zz_xy[i] * pa_x[i];

        ts_xzz_xz[i] = ts_zz_z[i] * fe_0 + ts_zz_xz[i] * pa_x[i];

        ts_xzz_yy[i] = ts_zz_yy[i] * pa_x[i];

        ts_xzz_yz[i] = ts_zz_yz[i] * pa_x[i];

        ts_xzz_zz[i] = ts_zz_zz[i] * pa_x[i];
    }

    // Set up 36-42 components of targeted buffer : FD

    auto ts_yyy_xx = pbuffer.data(idx_ovl_fd + 36);

    auto ts_yyy_xy = pbuffer.data(idx_ovl_fd + 37);

    auto ts_yyy_xz = pbuffer.data(idx_ovl_fd + 38);

    auto ts_yyy_yy = pbuffer.data(idx_ovl_fd + 39);

    auto ts_yyy_yz = pbuffer.data(idx_ovl_fd + 40);

    auto ts_yyy_zz = pbuffer.data(idx_ovl_fd + 41);

#pragma omp simd aligned(pa_y,          \
                             ts_y_xx,   \
                             ts_y_xy,   \
                             ts_y_xz,   \
                             ts_y_yy,   \
                             ts_y_yz,   \
                             ts_y_zz,   \
                             ts_yy_x,   \
                             ts_yy_xx,  \
                             ts_yy_xy,  \
                             ts_yy_xz,  \
                             ts_yy_y,   \
                             ts_yy_yy,  \
                             ts_yy_yz,  \
                             ts_yy_z,   \
                             ts_yy_zz,  \
                             ts_yyy_xx, \
                             ts_yyy_xy, \
                             ts_yyy_xz, \
                             ts_yyy_yy, \
                             ts_yyy_yz, \
                             ts_yyy_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyy_xx[i] = 2.0 * ts_y_xx[i] * fe_0 + ts_yy_xx[i] * pa_y[i];

        ts_yyy_xy[i] = 2.0 * ts_y_xy[i] * fe_0 + ts_yy_x[i] * fe_0 + ts_yy_xy[i] * pa_y[i];

        ts_yyy_xz[i] = 2.0 * ts_y_xz[i] * fe_0 + ts_yy_xz[i] * pa_y[i];

        ts_yyy_yy[i] = 2.0 * ts_y_yy[i] * fe_0 + 2.0 * ts_yy_y[i] * fe_0 + ts_yy_yy[i] * pa_y[i];

        ts_yyy_yz[i] = 2.0 * ts_y_yz[i] * fe_0 + ts_yy_z[i] * fe_0 + ts_yy_yz[i] * pa_y[i];

        ts_yyy_zz[i] = 2.0 * ts_y_zz[i] * fe_0 + ts_yy_zz[i] * pa_y[i];
    }

    // Set up 42-48 components of targeted buffer : FD

    auto ts_yyz_xx = pbuffer.data(idx_ovl_fd + 42);

    auto ts_yyz_xy = pbuffer.data(idx_ovl_fd + 43);

    auto ts_yyz_xz = pbuffer.data(idx_ovl_fd + 44);

    auto ts_yyz_yy = pbuffer.data(idx_ovl_fd + 45);

    auto ts_yyz_yz = pbuffer.data(idx_ovl_fd + 46);

    auto ts_yyz_zz = pbuffer.data(idx_ovl_fd + 47);

#pragma omp simd aligned(pa_y,          \
                             pa_z,      \
                             ts_yy_xx,  \
                             ts_yy_xy,  \
                             ts_yy_y,   \
                             ts_yy_yy,  \
                             ts_yy_yz,  \
                             ts_yyz_xx, \
                             ts_yyz_xy, \
                             ts_yyz_xz, \
                             ts_yyz_yy, \
                             ts_yyz_yz, \
                             ts_yyz_zz, \
                             ts_yz_xz,  \
                             ts_yz_zz,  \
                             ts_z_xz,   \
                             ts_z_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyz_xx[i] = ts_yy_xx[i] * pa_z[i];

        ts_yyz_xy[i] = ts_yy_xy[i] * pa_z[i];

        ts_yyz_xz[i] = ts_z_xz[i] * fe_0 + ts_yz_xz[i] * pa_y[i];

        ts_yyz_yy[i] = ts_yy_yy[i] * pa_z[i];

        ts_yyz_yz[i] = ts_yy_y[i] * fe_0 + ts_yy_yz[i] * pa_z[i];

        ts_yyz_zz[i] = ts_z_zz[i] * fe_0 + ts_yz_zz[i] * pa_y[i];
    }

    // Set up 48-54 components of targeted buffer : FD

    auto ts_yzz_xx = pbuffer.data(idx_ovl_fd + 48);

    auto ts_yzz_xy = pbuffer.data(idx_ovl_fd + 49);

    auto ts_yzz_xz = pbuffer.data(idx_ovl_fd + 50);

    auto ts_yzz_yy = pbuffer.data(idx_ovl_fd + 51);

    auto ts_yzz_yz = pbuffer.data(idx_ovl_fd + 52);

    auto ts_yzz_zz = pbuffer.data(idx_ovl_fd + 53);

#pragma omp simd aligned(pa_y,          \
                             ts_yzz_xx, \
                             ts_yzz_xy, \
                             ts_yzz_xz, \
                             ts_yzz_yy, \
                             ts_yzz_yz, \
                             ts_yzz_zz, \
                             ts_zz_x,   \
                             ts_zz_xx,  \
                             ts_zz_xy,  \
                             ts_zz_xz,  \
                             ts_zz_y,   \
                             ts_zz_yy,  \
                             ts_zz_yz,  \
                             ts_zz_z,   \
                             ts_zz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yzz_xx[i] = ts_zz_xx[i] * pa_y[i];

        ts_yzz_xy[i] = ts_zz_x[i] * fe_0 + ts_zz_xy[i] * pa_y[i];

        ts_yzz_xz[i] = ts_zz_xz[i] * pa_y[i];

        ts_yzz_yy[i] = 2.0 * ts_zz_y[i] * fe_0 + ts_zz_yy[i] * pa_y[i];

        ts_yzz_yz[i] = ts_zz_z[i] * fe_0 + ts_zz_yz[i] * pa_y[i];

        ts_yzz_zz[i] = ts_zz_zz[i] * pa_y[i];
    }

    // Set up 54-60 components of targeted buffer : FD

    auto ts_zzz_xx = pbuffer.data(idx_ovl_fd + 54);

    auto ts_zzz_xy = pbuffer.data(idx_ovl_fd + 55);

    auto ts_zzz_xz = pbuffer.data(idx_ovl_fd + 56);

    auto ts_zzz_yy = pbuffer.data(idx_ovl_fd + 57);

    auto ts_zzz_yz = pbuffer.data(idx_ovl_fd + 58);

    auto ts_zzz_zz = pbuffer.data(idx_ovl_fd + 59);

#pragma omp simd aligned(pa_z,          \
                             ts_z_xx,   \
                             ts_z_xy,   \
                             ts_z_xz,   \
                             ts_z_yy,   \
                             ts_z_yz,   \
                             ts_z_zz,   \
                             ts_zz_x,   \
                             ts_zz_xx,  \
                             ts_zz_xy,  \
                             ts_zz_xz,  \
                             ts_zz_y,   \
                             ts_zz_yy,  \
                             ts_zz_yz,  \
                             ts_zz_z,   \
                             ts_zz_zz,  \
                             ts_zzz_xx, \
                             ts_zzz_xy, \
                             ts_zzz_xz, \
                             ts_zzz_yy, \
                             ts_zzz_yz, \
                             ts_zzz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_zzz_xx[i] = 2.0 * ts_z_xx[i] * fe_0 + ts_zz_xx[i] * pa_z[i];

        ts_zzz_xy[i] = 2.0 * ts_z_xy[i] * fe_0 + ts_zz_xy[i] * pa_z[i];

        ts_zzz_xz[i] = 2.0 * ts_z_xz[i] * fe_0 + ts_zz_x[i] * fe_0 + ts_zz_xz[i] * pa_z[i];

        ts_zzz_yy[i] = 2.0 * ts_z_yy[i] * fe_0 + ts_zz_yy[i] * pa_z[i];

        ts_zzz_yz[i] = 2.0 * ts_z_yz[i] * fe_0 + ts_zz_y[i] * fe_0 + ts_zz_yz[i] * pa_z[i];

        ts_zzz_zz[i] = 2.0 * ts_z_zz[i] * fe_0 + 2.0 * ts_zz_z[i] * fe_0 + ts_zz_zz[i] * pa_z[i];
    }
}

}  // namespace ovlrec
