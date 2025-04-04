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

#include "OverlapPrimRecGD.hpp"

namespace ovlrec {  // ovlrec namespace

auto
comp_prim_overlap_gd(CSimdArray<double>&       pbuffer,
                     const size_t              idx_ovl_gd,
                     const size_t              idx_ovl_dd,
                     const size_t              idx_ovl_fp,
                     const size_t              idx_ovl_fd,
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

    // Set up components of auxiliary buffer : DD

    auto ts_xx_xx = pbuffer.data(idx_ovl_dd);

    auto ts_xx_xy = pbuffer.data(idx_ovl_dd + 1);

    auto ts_xx_xz = pbuffer.data(idx_ovl_dd + 2);

    auto ts_xx_yy = pbuffer.data(idx_ovl_dd + 3);

    auto ts_xx_yz = pbuffer.data(idx_ovl_dd + 4);

    auto ts_xx_zz = pbuffer.data(idx_ovl_dd + 5);

    auto ts_xy_yy = pbuffer.data(idx_ovl_dd + 9);

    auto ts_xy_yz = pbuffer.data(idx_ovl_dd + 10);

    auto ts_xz_yz = pbuffer.data(idx_ovl_dd + 16);

    auto ts_xz_zz = pbuffer.data(idx_ovl_dd + 17);

    auto ts_yy_xx = pbuffer.data(idx_ovl_dd + 18);

    auto ts_yy_xy = pbuffer.data(idx_ovl_dd + 19);

    auto ts_yy_xz = pbuffer.data(idx_ovl_dd + 20);

    auto ts_yy_yy = pbuffer.data(idx_ovl_dd + 21);

    auto ts_yy_yz = pbuffer.data(idx_ovl_dd + 22);

    auto ts_yy_zz = pbuffer.data(idx_ovl_dd + 23);

    auto ts_yz_xz = pbuffer.data(idx_ovl_dd + 26);

    auto ts_yz_yz = pbuffer.data(idx_ovl_dd + 28);

    auto ts_yz_zz = pbuffer.data(idx_ovl_dd + 29);

    auto ts_zz_xx = pbuffer.data(idx_ovl_dd + 30);

    auto ts_zz_xy = pbuffer.data(idx_ovl_dd + 31);

    auto ts_zz_xz = pbuffer.data(idx_ovl_dd + 32);

    auto ts_zz_yy = pbuffer.data(idx_ovl_dd + 33);

    auto ts_zz_yz = pbuffer.data(idx_ovl_dd + 34);

    auto ts_zz_zz = pbuffer.data(idx_ovl_dd + 35);

    // Set up components of auxiliary buffer : FP

    auto ts_xxx_x = pbuffer.data(idx_ovl_fp);

    auto ts_xxx_y = pbuffer.data(idx_ovl_fp + 1);

    auto ts_xxx_z = pbuffer.data(idx_ovl_fp + 2);

    auto ts_xyy_y = pbuffer.data(idx_ovl_fp + 10);

    auto ts_xzz_z = pbuffer.data(idx_ovl_fp + 17);

    auto ts_yyy_x = pbuffer.data(idx_ovl_fp + 18);

    auto ts_yyy_y = pbuffer.data(idx_ovl_fp + 19);

    auto ts_yyy_z = pbuffer.data(idx_ovl_fp + 20);

    auto ts_yyz_z = pbuffer.data(idx_ovl_fp + 23);

    auto ts_yzz_y = pbuffer.data(idx_ovl_fp + 25);

    auto ts_yzz_z = pbuffer.data(idx_ovl_fp + 26);

    auto ts_zzz_x = pbuffer.data(idx_ovl_fp + 27);

    auto ts_zzz_y = pbuffer.data(idx_ovl_fp + 28);

    auto ts_zzz_z = pbuffer.data(idx_ovl_fp + 29);

    // Set up components of auxiliary buffer : FD

    auto ts_xxx_xx = pbuffer.data(idx_ovl_fd);

    auto ts_xxx_xy = pbuffer.data(idx_ovl_fd + 1);

    auto ts_xxx_xz = pbuffer.data(idx_ovl_fd + 2);

    auto ts_xxx_yy = pbuffer.data(idx_ovl_fd + 3);

    auto ts_xxx_yz = pbuffer.data(idx_ovl_fd + 4);

    auto ts_xxx_zz = pbuffer.data(idx_ovl_fd + 5);

    auto ts_xxy_xx = pbuffer.data(idx_ovl_fd + 6);

    auto ts_xxy_xy = pbuffer.data(idx_ovl_fd + 7);

    auto ts_xxy_xz = pbuffer.data(idx_ovl_fd + 8);

    auto ts_xxy_yy = pbuffer.data(idx_ovl_fd + 9);

    auto ts_xxy_yz = pbuffer.data(idx_ovl_fd + 10);

    auto ts_xxz_xx = pbuffer.data(idx_ovl_fd + 12);

    auto ts_xxz_xy = pbuffer.data(idx_ovl_fd + 13);

    auto ts_xxz_xz = pbuffer.data(idx_ovl_fd + 14);

    auto ts_xxz_yz = pbuffer.data(idx_ovl_fd + 16);

    auto ts_xxz_zz = pbuffer.data(idx_ovl_fd + 17);

    auto ts_xyy_xx = pbuffer.data(idx_ovl_fd + 18);

    auto ts_xyy_xy = pbuffer.data(idx_ovl_fd + 19);

    auto ts_xyy_yy = pbuffer.data(idx_ovl_fd + 21);

    auto ts_xyy_yz = pbuffer.data(idx_ovl_fd + 22);

    auto ts_xyy_zz = pbuffer.data(idx_ovl_fd + 23);

    auto ts_xyz_yz = pbuffer.data(idx_ovl_fd + 28);

    auto ts_xzz_xx = pbuffer.data(idx_ovl_fd + 30);

    auto ts_xzz_xz = pbuffer.data(idx_ovl_fd + 32);

    auto ts_xzz_yy = pbuffer.data(idx_ovl_fd + 33);

    auto ts_xzz_yz = pbuffer.data(idx_ovl_fd + 34);

    auto ts_xzz_zz = pbuffer.data(idx_ovl_fd + 35);

    auto ts_yyy_xx = pbuffer.data(idx_ovl_fd + 36);

    auto ts_yyy_xy = pbuffer.data(idx_ovl_fd + 37);

    auto ts_yyy_xz = pbuffer.data(idx_ovl_fd + 38);

    auto ts_yyy_yy = pbuffer.data(idx_ovl_fd + 39);

    auto ts_yyy_yz = pbuffer.data(idx_ovl_fd + 40);

    auto ts_yyy_zz = pbuffer.data(idx_ovl_fd + 41);

    auto ts_yyz_xy = pbuffer.data(idx_ovl_fd + 43);

    auto ts_yyz_xz = pbuffer.data(idx_ovl_fd + 44);

    auto ts_yyz_yy = pbuffer.data(idx_ovl_fd + 45);

    auto ts_yyz_yz = pbuffer.data(idx_ovl_fd + 46);

    auto ts_yyz_zz = pbuffer.data(idx_ovl_fd + 47);

    auto ts_yzz_xx = pbuffer.data(idx_ovl_fd + 48);

    auto ts_yzz_xy = pbuffer.data(idx_ovl_fd + 49);

    auto ts_yzz_xz = pbuffer.data(idx_ovl_fd + 50);

    auto ts_yzz_yy = pbuffer.data(idx_ovl_fd + 51);

    auto ts_yzz_yz = pbuffer.data(idx_ovl_fd + 52);

    auto ts_yzz_zz = pbuffer.data(idx_ovl_fd + 53);

    auto ts_zzz_xx = pbuffer.data(idx_ovl_fd + 54);

    auto ts_zzz_xy = pbuffer.data(idx_ovl_fd + 55);

    auto ts_zzz_xz = pbuffer.data(idx_ovl_fd + 56);

    auto ts_zzz_yy = pbuffer.data(idx_ovl_fd + 57);

    auto ts_zzz_yz = pbuffer.data(idx_ovl_fd + 58);

    auto ts_zzz_zz = pbuffer.data(idx_ovl_fd + 59);

    // Set up 0-6 components of targeted buffer : GD

    auto ts_xxxx_xx = pbuffer.data(idx_ovl_gd);

    auto ts_xxxx_xy = pbuffer.data(idx_ovl_gd + 1);

    auto ts_xxxx_xz = pbuffer.data(idx_ovl_gd + 2);

    auto ts_xxxx_yy = pbuffer.data(idx_ovl_gd + 3);

    auto ts_xxxx_yz = pbuffer.data(idx_ovl_gd + 4);

    auto ts_xxxx_zz = pbuffer.data(idx_ovl_gd + 5);

#pragma omp simd aligned(pa_x,           \
                             ts_xx_xx,   \
                             ts_xx_xy,   \
                             ts_xx_xz,   \
                             ts_xx_yy,   \
                             ts_xx_yz,   \
                             ts_xx_zz,   \
                             ts_xxx_x,   \
                             ts_xxx_xx,  \
                             ts_xxx_xy,  \
                             ts_xxx_xz,  \
                             ts_xxx_y,   \
                             ts_xxx_yy,  \
                             ts_xxx_yz,  \
                             ts_xxx_z,   \
                             ts_xxx_zz,  \
                             ts_xxxx_xx, \
                             ts_xxxx_xy, \
                             ts_xxxx_xz, \
                             ts_xxxx_yy, \
                             ts_xxxx_yz, \
                             ts_xxxx_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxx_xx[i] = 3.0 * ts_xx_xx[i] * fe_0 + 2.0 * ts_xxx_x[i] * fe_0 + ts_xxx_xx[i] * pa_x[i];

        ts_xxxx_xy[i] = 3.0 * ts_xx_xy[i] * fe_0 + ts_xxx_y[i] * fe_0 + ts_xxx_xy[i] * pa_x[i];

        ts_xxxx_xz[i] = 3.0 * ts_xx_xz[i] * fe_0 + ts_xxx_z[i] * fe_0 + ts_xxx_xz[i] * pa_x[i];

        ts_xxxx_yy[i] = 3.0 * ts_xx_yy[i] * fe_0 + ts_xxx_yy[i] * pa_x[i];

        ts_xxxx_yz[i] = 3.0 * ts_xx_yz[i] * fe_0 + ts_xxx_yz[i] * pa_x[i];

        ts_xxxx_zz[i] = 3.0 * ts_xx_zz[i] * fe_0 + ts_xxx_zz[i] * pa_x[i];
    }

    // Set up 6-12 components of targeted buffer : GD

    auto ts_xxxy_xx = pbuffer.data(idx_ovl_gd + 6);

    auto ts_xxxy_xy = pbuffer.data(idx_ovl_gd + 7);

    auto ts_xxxy_xz = pbuffer.data(idx_ovl_gd + 8);

    auto ts_xxxy_yy = pbuffer.data(idx_ovl_gd + 9);

    auto ts_xxxy_yz = pbuffer.data(idx_ovl_gd + 10);

    auto ts_xxxy_zz = pbuffer.data(idx_ovl_gd + 11);

#pragma omp simd aligned(pa_x,           \
                             pa_y,       \
                             ts_xxx_x,   \
                             ts_xxx_xx,  \
                             ts_xxx_xy,  \
                             ts_xxx_xz,  \
                             ts_xxx_zz,  \
                             ts_xxxy_xx, \
                             ts_xxxy_xy, \
                             ts_xxxy_xz, \
                             ts_xxxy_yy, \
                             ts_xxxy_yz, \
                             ts_xxxy_zz, \
                             ts_xxy_yy,  \
                             ts_xxy_yz,  \
                             ts_xy_yy,   \
                             ts_xy_yz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxy_xx[i] = ts_xxx_xx[i] * pa_y[i];

        ts_xxxy_xy[i] = ts_xxx_x[i] * fe_0 + ts_xxx_xy[i] * pa_y[i];

        ts_xxxy_xz[i] = ts_xxx_xz[i] * pa_y[i];

        ts_xxxy_yy[i] = 2.0 * ts_xy_yy[i] * fe_0 + ts_xxy_yy[i] * pa_x[i];

        ts_xxxy_yz[i] = 2.0 * ts_xy_yz[i] * fe_0 + ts_xxy_yz[i] * pa_x[i];

        ts_xxxy_zz[i] = ts_xxx_zz[i] * pa_y[i];
    }

    // Set up 12-18 components of targeted buffer : GD

    auto ts_xxxz_xx = pbuffer.data(idx_ovl_gd + 12);

    auto ts_xxxz_xy = pbuffer.data(idx_ovl_gd + 13);

    auto ts_xxxz_xz = pbuffer.data(idx_ovl_gd + 14);

    auto ts_xxxz_yy = pbuffer.data(idx_ovl_gd + 15);

    auto ts_xxxz_yz = pbuffer.data(idx_ovl_gd + 16);

    auto ts_xxxz_zz = pbuffer.data(idx_ovl_gd + 17);

#pragma omp simd aligned(pa_x,           \
                             pa_z,       \
                             ts_xxx_x,   \
                             ts_xxx_xx,  \
                             ts_xxx_xy,  \
                             ts_xxx_xz,  \
                             ts_xxx_yy,  \
                             ts_xxxz_xx, \
                             ts_xxxz_xy, \
                             ts_xxxz_xz, \
                             ts_xxxz_yy, \
                             ts_xxxz_yz, \
                             ts_xxxz_zz, \
                             ts_xxz_yz,  \
                             ts_xxz_zz,  \
                             ts_xz_yz,   \
                             ts_xz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxz_xx[i] = ts_xxx_xx[i] * pa_z[i];

        ts_xxxz_xy[i] = ts_xxx_xy[i] * pa_z[i];

        ts_xxxz_xz[i] = ts_xxx_x[i] * fe_0 + ts_xxx_xz[i] * pa_z[i];

        ts_xxxz_yy[i] = ts_xxx_yy[i] * pa_z[i];

        ts_xxxz_yz[i] = 2.0 * ts_xz_yz[i] * fe_0 + ts_xxz_yz[i] * pa_x[i];

        ts_xxxz_zz[i] = 2.0 * ts_xz_zz[i] * fe_0 + ts_xxz_zz[i] * pa_x[i];
    }

    // Set up 18-24 components of targeted buffer : GD

    auto ts_xxyy_xx = pbuffer.data(idx_ovl_gd + 18);

    auto ts_xxyy_xy = pbuffer.data(idx_ovl_gd + 19);

    auto ts_xxyy_xz = pbuffer.data(idx_ovl_gd + 20);

    auto ts_xxyy_yy = pbuffer.data(idx_ovl_gd + 21);

    auto ts_xxyy_yz = pbuffer.data(idx_ovl_gd + 22);

    auto ts_xxyy_zz = pbuffer.data(idx_ovl_gd + 23);

#pragma omp simd aligned(pa_x,           \
                             pa_y,       \
                             ts_xx_xx,   \
                             ts_xx_xz,   \
                             ts_xxy_xx,  \
                             ts_xxy_xz,  \
                             ts_xxyy_xx, \
                             ts_xxyy_xy, \
                             ts_xxyy_xz, \
                             ts_xxyy_yy, \
                             ts_xxyy_yz, \
                             ts_xxyy_zz, \
                             ts_xyy_xy,  \
                             ts_xyy_y,   \
                             ts_xyy_yy,  \
                             ts_xyy_yz,  \
                             ts_xyy_zz,  \
                             ts_yy_xy,   \
                             ts_yy_yy,   \
                             ts_yy_yz,   \
                             ts_yy_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyy_xx[i] = ts_xx_xx[i] * fe_0 + ts_xxy_xx[i] * pa_y[i];

        ts_xxyy_xy[i] = ts_yy_xy[i] * fe_0 + ts_xyy_y[i] * fe_0 + ts_xyy_xy[i] * pa_x[i];

        ts_xxyy_xz[i] = ts_xx_xz[i] * fe_0 + ts_xxy_xz[i] * pa_y[i];

        ts_xxyy_yy[i] = ts_yy_yy[i] * fe_0 + ts_xyy_yy[i] * pa_x[i];

        ts_xxyy_yz[i] = ts_yy_yz[i] * fe_0 + ts_xyy_yz[i] * pa_x[i];

        ts_xxyy_zz[i] = ts_yy_zz[i] * fe_0 + ts_xyy_zz[i] * pa_x[i];
    }

    // Set up 24-30 components of targeted buffer : GD

    auto ts_xxyz_xx = pbuffer.data(idx_ovl_gd + 24);

    auto ts_xxyz_xy = pbuffer.data(idx_ovl_gd + 25);

    auto ts_xxyz_xz = pbuffer.data(idx_ovl_gd + 26);

    auto ts_xxyz_yy = pbuffer.data(idx_ovl_gd + 27);

    auto ts_xxyz_yz = pbuffer.data(idx_ovl_gd + 28);

    auto ts_xxyz_zz = pbuffer.data(idx_ovl_gd + 29);

#pragma omp simd aligned(pa_x,           \
                             pa_y,       \
                             pa_z,       \
                             ts_xxy_xy,  \
                             ts_xxy_yy,  \
                             ts_xxyz_xx, \
                             ts_xxyz_xy, \
                             ts_xxyz_xz, \
                             ts_xxyz_yy, \
                             ts_xxyz_yz, \
                             ts_xxyz_zz, \
                             ts_xxz_xx,  \
                             ts_xxz_xz,  \
                             ts_xxz_zz,  \
                             ts_xyz_yz,  \
                             ts_yz_yz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyz_xx[i] = ts_xxz_xx[i] * pa_y[i];

        ts_xxyz_xy[i] = ts_xxy_xy[i] * pa_z[i];

        ts_xxyz_xz[i] = ts_xxz_xz[i] * pa_y[i];

        ts_xxyz_yy[i] = ts_xxy_yy[i] * pa_z[i];

        ts_xxyz_yz[i] = ts_yz_yz[i] * fe_0 + ts_xyz_yz[i] * pa_x[i];

        ts_xxyz_zz[i] = ts_xxz_zz[i] * pa_y[i];
    }

    // Set up 30-36 components of targeted buffer : GD

    auto ts_xxzz_xx = pbuffer.data(idx_ovl_gd + 30);

    auto ts_xxzz_xy = pbuffer.data(idx_ovl_gd + 31);

    auto ts_xxzz_xz = pbuffer.data(idx_ovl_gd + 32);

    auto ts_xxzz_yy = pbuffer.data(idx_ovl_gd + 33);

    auto ts_xxzz_yz = pbuffer.data(idx_ovl_gd + 34);

    auto ts_xxzz_zz = pbuffer.data(idx_ovl_gd + 35);

#pragma omp simd aligned(pa_x,           \
                             pa_z,       \
                             ts_xx_xx,   \
                             ts_xx_xy,   \
                             ts_xxz_xx,  \
                             ts_xxz_xy,  \
                             ts_xxzz_xx, \
                             ts_xxzz_xy, \
                             ts_xxzz_xz, \
                             ts_xxzz_yy, \
                             ts_xxzz_yz, \
                             ts_xxzz_zz, \
                             ts_xzz_xz,  \
                             ts_xzz_yy,  \
                             ts_xzz_yz,  \
                             ts_xzz_z,   \
                             ts_xzz_zz,  \
                             ts_zz_xz,   \
                             ts_zz_yy,   \
                             ts_zz_yz,   \
                             ts_zz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxzz_xx[i] = ts_xx_xx[i] * fe_0 + ts_xxz_xx[i] * pa_z[i];

        ts_xxzz_xy[i] = ts_xx_xy[i] * fe_0 + ts_xxz_xy[i] * pa_z[i];

        ts_xxzz_xz[i] = ts_zz_xz[i] * fe_0 + ts_xzz_z[i] * fe_0 + ts_xzz_xz[i] * pa_x[i];

        ts_xxzz_yy[i] = ts_zz_yy[i] * fe_0 + ts_xzz_yy[i] * pa_x[i];

        ts_xxzz_yz[i] = ts_zz_yz[i] * fe_0 + ts_xzz_yz[i] * pa_x[i];

        ts_xxzz_zz[i] = ts_zz_zz[i] * fe_0 + ts_xzz_zz[i] * pa_x[i];
    }

    // Set up 36-42 components of targeted buffer : GD

    auto ts_xyyy_xx = pbuffer.data(idx_ovl_gd + 36);

    auto ts_xyyy_xy = pbuffer.data(idx_ovl_gd + 37);

    auto ts_xyyy_xz = pbuffer.data(idx_ovl_gd + 38);

    auto ts_xyyy_yy = pbuffer.data(idx_ovl_gd + 39);

    auto ts_xyyy_yz = pbuffer.data(idx_ovl_gd + 40);

    auto ts_xyyy_zz = pbuffer.data(idx_ovl_gd + 41);

#pragma omp simd aligned(pa_x,           \
                             ts_xyyy_xx, \
                             ts_xyyy_xy, \
                             ts_xyyy_xz, \
                             ts_xyyy_yy, \
                             ts_xyyy_yz, \
                             ts_xyyy_zz, \
                             ts_yyy_x,   \
                             ts_yyy_xx,  \
                             ts_yyy_xy,  \
                             ts_yyy_xz,  \
                             ts_yyy_y,   \
                             ts_yyy_yy,  \
                             ts_yyy_yz,  \
                             ts_yyy_z,   \
                             ts_yyy_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyy_xx[i] = 2.0 * ts_yyy_x[i] * fe_0 + ts_yyy_xx[i] * pa_x[i];

        ts_xyyy_xy[i] = ts_yyy_y[i] * fe_0 + ts_yyy_xy[i] * pa_x[i];

        ts_xyyy_xz[i] = ts_yyy_z[i] * fe_0 + ts_yyy_xz[i] * pa_x[i];

        ts_xyyy_yy[i] = ts_yyy_yy[i] * pa_x[i];

        ts_xyyy_yz[i] = ts_yyy_yz[i] * pa_x[i];

        ts_xyyy_zz[i] = ts_yyy_zz[i] * pa_x[i];
    }

    // Set up 42-48 components of targeted buffer : GD

    auto ts_xyyz_xx = pbuffer.data(idx_ovl_gd + 42);

    auto ts_xyyz_xy = pbuffer.data(idx_ovl_gd + 43);

    auto ts_xyyz_xz = pbuffer.data(idx_ovl_gd + 44);

    auto ts_xyyz_yy = pbuffer.data(idx_ovl_gd + 45);

    auto ts_xyyz_yz = pbuffer.data(idx_ovl_gd + 46);

    auto ts_xyyz_zz = pbuffer.data(idx_ovl_gd + 47);

#pragma omp simd aligned(pa_x,           \
                             pa_z,       \
                             ts_xyy_xx,  \
                             ts_xyy_xy,  \
                             ts_xyyz_xx, \
                             ts_xyyz_xy, \
                             ts_xyyz_xz, \
                             ts_xyyz_yy, \
                             ts_xyyz_yz, \
                             ts_xyyz_zz, \
                             ts_yyz_xz,  \
                             ts_yyz_yy,  \
                             ts_yyz_yz,  \
                             ts_yyz_z,   \
                             ts_yyz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyz_xx[i] = ts_xyy_xx[i] * pa_z[i];

        ts_xyyz_xy[i] = ts_xyy_xy[i] * pa_z[i];

        ts_xyyz_xz[i] = ts_yyz_z[i] * fe_0 + ts_yyz_xz[i] * pa_x[i];

        ts_xyyz_yy[i] = ts_yyz_yy[i] * pa_x[i];

        ts_xyyz_yz[i] = ts_yyz_yz[i] * pa_x[i];

        ts_xyyz_zz[i] = ts_yyz_zz[i] * pa_x[i];
    }

    // Set up 48-54 components of targeted buffer : GD

    auto ts_xyzz_xx = pbuffer.data(idx_ovl_gd + 48);

    auto ts_xyzz_xy = pbuffer.data(idx_ovl_gd + 49);

    auto ts_xyzz_xz = pbuffer.data(idx_ovl_gd + 50);

    auto ts_xyzz_yy = pbuffer.data(idx_ovl_gd + 51);

    auto ts_xyzz_yz = pbuffer.data(idx_ovl_gd + 52);

    auto ts_xyzz_zz = pbuffer.data(idx_ovl_gd + 53);

#pragma omp simd aligned(pa_x,           \
                             pa_y,       \
                             ts_xyzz_xx, \
                             ts_xyzz_xy, \
                             ts_xyzz_xz, \
                             ts_xyzz_yy, \
                             ts_xyzz_yz, \
                             ts_xyzz_zz, \
                             ts_xzz_xx,  \
                             ts_xzz_xz,  \
                             ts_yzz_xy,  \
                             ts_yzz_y,   \
                             ts_yzz_yy,  \
                             ts_yzz_yz,  \
                             ts_yzz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyzz_xx[i] = ts_xzz_xx[i] * pa_y[i];

        ts_xyzz_xy[i] = ts_yzz_y[i] * fe_0 + ts_yzz_xy[i] * pa_x[i];

        ts_xyzz_xz[i] = ts_xzz_xz[i] * pa_y[i];

        ts_xyzz_yy[i] = ts_yzz_yy[i] * pa_x[i];

        ts_xyzz_yz[i] = ts_yzz_yz[i] * pa_x[i];

        ts_xyzz_zz[i] = ts_yzz_zz[i] * pa_x[i];
    }

    // Set up 54-60 components of targeted buffer : GD

    auto ts_xzzz_xx = pbuffer.data(idx_ovl_gd + 54);

    auto ts_xzzz_xy = pbuffer.data(idx_ovl_gd + 55);

    auto ts_xzzz_xz = pbuffer.data(idx_ovl_gd + 56);

    auto ts_xzzz_yy = pbuffer.data(idx_ovl_gd + 57);

    auto ts_xzzz_yz = pbuffer.data(idx_ovl_gd + 58);

    auto ts_xzzz_zz = pbuffer.data(idx_ovl_gd + 59);

#pragma omp simd aligned(pa_x,           \
                             ts_xzzz_xx, \
                             ts_xzzz_xy, \
                             ts_xzzz_xz, \
                             ts_xzzz_yy, \
                             ts_xzzz_yz, \
                             ts_xzzz_zz, \
                             ts_zzz_x,   \
                             ts_zzz_xx,  \
                             ts_zzz_xy,  \
                             ts_zzz_xz,  \
                             ts_zzz_y,   \
                             ts_zzz_yy,  \
                             ts_zzz_yz,  \
                             ts_zzz_z,   \
                             ts_zzz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xzzz_xx[i] = 2.0 * ts_zzz_x[i] * fe_0 + ts_zzz_xx[i] * pa_x[i];

        ts_xzzz_xy[i] = ts_zzz_y[i] * fe_0 + ts_zzz_xy[i] * pa_x[i];

        ts_xzzz_xz[i] = ts_zzz_z[i] * fe_0 + ts_zzz_xz[i] * pa_x[i];

        ts_xzzz_yy[i] = ts_zzz_yy[i] * pa_x[i];

        ts_xzzz_yz[i] = ts_zzz_yz[i] * pa_x[i];

        ts_xzzz_zz[i] = ts_zzz_zz[i] * pa_x[i];
    }

    // Set up 60-66 components of targeted buffer : GD

    auto ts_yyyy_xx = pbuffer.data(idx_ovl_gd + 60);

    auto ts_yyyy_xy = pbuffer.data(idx_ovl_gd + 61);

    auto ts_yyyy_xz = pbuffer.data(idx_ovl_gd + 62);

    auto ts_yyyy_yy = pbuffer.data(idx_ovl_gd + 63);

    auto ts_yyyy_yz = pbuffer.data(idx_ovl_gd + 64);

    auto ts_yyyy_zz = pbuffer.data(idx_ovl_gd + 65);

#pragma omp simd aligned(pa_y,           \
                             ts_yy_xx,   \
                             ts_yy_xy,   \
                             ts_yy_xz,   \
                             ts_yy_yy,   \
                             ts_yy_yz,   \
                             ts_yy_zz,   \
                             ts_yyy_x,   \
                             ts_yyy_xx,  \
                             ts_yyy_xy,  \
                             ts_yyy_xz,  \
                             ts_yyy_y,   \
                             ts_yyy_yy,  \
                             ts_yyy_yz,  \
                             ts_yyy_z,   \
                             ts_yyy_zz,  \
                             ts_yyyy_xx, \
                             ts_yyyy_xy, \
                             ts_yyyy_xz, \
                             ts_yyyy_yy, \
                             ts_yyyy_yz, \
                             ts_yyyy_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyy_xx[i] = 3.0 * ts_yy_xx[i] * fe_0 + ts_yyy_xx[i] * pa_y[i];

        ts_yyyy_xy[i] = 3.0 * ts_yy_xy[i] * fe_0 + ts_yyy_x[i] * fe_0 + ts_yyy_xy[i] * pa_y[i];

        ts_yyyy_xz[i] = 3.0 * ts_yy_xz[i] * fe_0 + ts_yyy_xz[i] * pa_y[i];

        ts_yyyy_yy[i] = 3.0 * ts_yy_yy[i] * fe_0 + 2.0 * ts_yyy_y[i] * fe_0 + ts_yyy_yy[i] * pa_y[i];

        ts_yyyy_yz[i] = 3.0 * ts_yy_yz[i] * fe_0 + ts_yyy_z[i] * fe_0 + ts_yyy_yz[i] * pa_y[i];

        ts_yyyy_zz[i] = 3.0 * ts_yy_zz[i] * fe_0 + ts_yyy_zz[i] * pa_y[i];
    }

    // Set up 66-72 components of targeted buffer : GD

    auto ts_yyyz_xx = pbuffer.data(idx_ovl_gd + 66);

    auto ts_yyyz_xy = pbuffer.data(idx_ovl_gd + 67);

    auto ts_yyyz_xz = pbuffer.data(idx_ovl_gd + 68);

    auto ts_yyyz_yy = pbuffer.data(idx_ovl_gd + 69);

    auto ts_yyyz_yz = pbuffer.data(idx_ovl_gd + 70);

    auto ts_yyyz_zz = pbuffer.data(idx_ovl_gd + 71);

#pragma omp simd aligned(pa_y,           \
                             pa_z,       \
                             ts_yyy_xx,  \
                             ts_yyy_xy,  \
                             ts_yyy_y,   \
                             ts_yyy_yy,  \
                             ts_yyy_yz,  \
                             ts_yyyz_xx, \
                             ts_yyyz_xy, \
                             ts_yyyz_xz, \
                             ts_yyyz_yy, \
                             ts_yyyz_yz, \
                             ts_yyyz_zz, \
                             ts_yyz_xz,  \
                             ts_yyz_zz,  \
                             ts_yz_xz,   \
                             ts_yz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyz_xx[i] = ts_yyy_xx[i] * pa_z[i];

        ts_yyyz_xy[i] = ts_yyy_xy[i] * pa_z[i];

        ts_yyyz_xz[i] = 2.0 * ts_yz_xz[i] * fe_0 + ts_yyz_xz[i] * pa_y[i];

        ts_yyyz_yy[i] = ts_yyy_yy[i] * pa_z[i];

        ts_yyyz_yz[i] = ts_yyy_y[i] * fe_0 + ts_yyy_yz[i] * pa_z[i];

        ts_yyyz_zz[i] = 2.0 * ts_yz_zz[i] * fe_0 + ts_yyz_zz[i] * pa_y[i];
    }

    // Set up 72-78 components of targeted buffer : GD

    auto ts_yyzz_xx = pbuffer.data(idx_ovl_gd + 72);

    auto ts_yyzz_xy = pbuffer.data(idx_ovl_gd + 73);

    auto ts_yyzz_xz = pbuffer.data(idx_ovl_gd + 74);

    auto ts_yyzz_yy = pbuffer.data(idx_ovl_gd + 75);

    auto ts_yyzz_yz = pbuffer.data(idx_ovl_gd + 76);

    auto ts_yyzz_zz = pbuffer.data(idx_ovl_gd + 77);

#pragma omp simd aligned(pa_y,           \
                             pa_z,       \
                             ts_yy_xy,   \
                             ts_yy_yy,   \
                             ts_yyz_xy,  \
                             ts_yyz_yy,  \
                             ts_yyzz_xx, \
                             ts_yyzz_xy, \
                             ts_yyzz_xz, \
                             ts_yyzz_yy, \
                             ts_yyzz_yz, \
                             ts_yyzz_zz, \
                             ts_yzz_xx,  \
                             ts_yzz_xz,  \
                             ts_yzz_yz,  \
                             ts_yzz_z,   \
                             ts_yzz_zz,  \
                             ts_zz_xx,   \
                             ts_zz_xz,   \
                             ts_zz_yz,   \
                             ts_zz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyzz_xx[i] = ts_zz_xx[i] * fe_0 + ts_yzz_xx[i] * pa_y[i];

        ts_yyzz_xy[i] = ts_yy_xy[i] * fe_0 + ts_yyz_xy[i] * pa_z[i];

        ts_yyzz_xz[i] = ts_zz_xz[i] * fe_0 + ts_yzz_xz[i] * pa_y[i];

        ts_yyzz_yy[i] = ts_yy_yy[i] * fe_0 + ts_yyz_yy[i] * pa_z[i];

        ts_yyzz_yz[i] = ts_zz_yz[i] * fe_0 + ts_yzz_z[i] * fe_0 + ts_yzz_yz[i] * pa_y[i];

        ts_yyzz_zz[i] = ts_zz_zz[i] * fe_0 + ts_yzz_zz[i] * pa_y[i];
    }

    // Set up 78-84 components of targeted buffer : GD

    auto ts_yzzz_xx = pbuffer.data(idx_ovl_gd + 78);

    auto ts_yzzz_xy = pbuffer.data(idx_ovl_gd + 79);

    auto ts_yzzz_xz = pbuffer.data(idx_ovl_gd + 80);

    auto ts_yzzz_yy = pbuffer.data(idx_ovl_gd + 81);

    auto ts_yzzz_yz = pbuffer.data(idx_ovl_gd + 82);

    auto ts_yzzz_zz = pbuffer.data(idx_ovl_gd + 83);

#pragma omp simd aligned(pa_y,           \
                             ts_yzzz_xx, \
                             ts_yzzz_xy, \
                             ts_yzzz_xz, \
                             ts_yzzz_yy, \
                             ts_yzzz_yz, \
                             ts_yzzz_zz, \
                             ts_zzz_x,   \
                             ts_zzz_xx,  \
                             ts_zzz_xy,  \
                             ts_zzz_xz,  \
                             ts_zzz_y,   \
                             ts_zzz_yy,  \
                             ts_zzz_yz,  \
                             ts_zzz_z,   \
                             ts_zzz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yzzz_xx[i] = ts_zzz_xx[i] * pa_y[i];

        ts_yzzz_xy[i] = ts_zzz_x[i] * fe_0 + ts_zzz_xy[i] * pa_y[i];

        ts_yzzz_xz[i] = ts_zzz_xz[i] * pa_y[i];

        ts_yzzz_yy[i] = 2.0 * ts_zzz_y[i] * fe_0 + ts_zzz_yy[i] * pa_y[i];

        ts_yzzz_yz[i] = ts_zzz_z[i] * fe_0 + ts_zzz_yz[i] * pa_y[i];

        ts_yzzz_zz[i] = ts_zzz_zz[i] * pa_y[i];
    }

    // Set up 84-90 components of targeted buffer : GD

    auto ts_zzzz_xx = pbuffer.data(idx_ovl_gd + 84);

    auto ts_zzzz_xy = pbuffer.data(idx_ovl_gd + 85);

    auto ts_zzzz_xz = pbuffer.data(idx_ovl_gd + 86);

    auto ts_zzzz_yy = pbuffer.data(idx_ovl_gd + 87);

    auto ts_zzzz_yz = pbuffer.data(idx_ovl_gd + 88);

    auto ts_zzzz_zz = pbuffer.data(idx_ovl_gd + 89);

#pragma omp simd aligned(pa_z,           \
                             ts_zz_xx,   \
                             ts_zz_xy,   \
                             ts_zz_xz,   \
                             ts_zz_yy,   \
                             ts_zz_yz,   \
                             ts_zz_zz,   \
                             ts_zzz_x,   \
                             ts_zzz_xx,  \
                             ts_zzz_xy,  \
                             ts_zzz_xz,  \
                             ts_zzz_y,   \
                             ts_zzz_yy,  \
                             ts_zzz_yz,  \
                             ts_zzz_z,   \
                             ts_zzz_zz,  \
                             ts_zzzz_xx, \
                             ts_zzzz_xy, \
                             ts_zzzz_xz, \
                             ts_zzzz_yy, \
                             ts_zzzz_yz, \
                             ts_zzzz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_zzzz_xx[i] = 3.0 * ts_zz_xx[i] * fe_0 + ts_zzz_xx[i] * pa_z[i];

        ts_zzzz_xy[i] = 3.0 * ts_zz_xy[i] * fe_0 + ts_zzz_xy[i] * pa_z[i];

        ts_zzzz_xz[i] = 3.0 * ts_zz_xz[i] * fe_0 + ts_zzz_x[i] * fe_0 + ts_zzz_xz[i] * pa_z[i];

        ts_zzzz_yy[i] = 3.0 * ts_zz_yy[i] * fe_0 + ts_zzz_yy[i] * pa_z[i];

        ts_zzzz_yz[i] = 3.0 * ts_zz_yz[i] * fe_0 + ts_zzz_y[i] * fe_0 + ts_zzz_yz[i] * pa_z[i];

        ts_zzzz_zz[i] = 3.0 * ts_zz_zz[i] * fe_0 + 2.0 * ts_zzz_z[i] * fe_0 + ts_zzz_zz[i] * pa_z[i];
    }
}

}  // namespace ovlrec
