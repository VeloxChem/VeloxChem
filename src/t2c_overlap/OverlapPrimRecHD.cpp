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

#include "OverlapPrimRecHD.hpp"

namespace ovlrec { // ovlrec namespace

auto
comp_prim_overlap_hd(CSimdArray<double>& pbuffer, 
                     const size_t idx_ovl_hd,
                     const size_t idx_ovl_fd,
                     const size_t idx_ovl_gp,
                     const size_t idx_ovl_gd,
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

    // Set up components of auxiliary buffer : FD

    auto ts_xxx_xx = pbuffer.data(idx_ovl_fd);

    auto ts_xxx_xy = pbuffer.data(idx_ovl_fd + 1);

    auto ts_xxx_xz = pbuffer.data(idx_ovl_fd + 2);

    auto ts_xxx_yy = pbuffer.data(idx_ovl_fd + 3);

    auto ts_xxx_yz = pbuffer.data(idx_ovl_fd + 4);

    auto ts_xxx_zz = pbuffer.data(idx_ovl_fd + 5);

    auto ts_xxy_xx = pbuffer.data(idx_ovl_fd + 6);

    auto ts_xxy_xz = pbuffer.data(idx_ovl_fd + 8);

    auto ts_xxy_yy = pbuffer.data(idx_ovl_fd + 9);

    auto ts_xxy_yz = pbuffer.data(idx_ovl_fd + 10);

    auto ts_xxz_xx = pbuffer.data(idx_ovl_fd + 12);

    auto ts_xxz_xy = pbuffer.data(idx_ovl_fd + 13);

    auto ts_xxz_xz = pbuffer.data(idx_ovl_fd + 14);

    auto ts_xxz_yz = pbuffer.data(idx_ovl_fd + 16);

    auto ts_xxz_zz = pbuffer.data(idx_ovl_fd + 17);

    auto ts_xyy_xy = pbuffer.data(idx_ovl_fd + 19);

    auto ts_xyy_yy = pbuffer.data(idx_ovl_fd + 21);

    auto ts_xyy_yz = pbuffer.data(idx_ovl_fd + 22);

    auto ts_xyy_zz = pbuffer.data(idx_ovl_fd + 23);

    auto ts_xyz_yz = pbuffer.data(idx_ovl_fd + 28);

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

    // Set up components of auxiliary buffer : GP

    auto ts_xxxx_x = pbuffer.data(idx_ovl_gp);

    auto ts_xxxx_y = pbuffer.data(idx_ovl_gp + 1);

    auto ts_xxxx_z = pbuffer.data(idx_ovl_gp + 2);

    auto ts_xxyy_y = pbuffer.data(idx_ovl_gp + 10);

    auto ts_xxzz_x = pbuffer.data(idx_ovl_gp + 15);

    auto ts_xxzz_z = pbuffer.data(idx_ovl_gp + 17);

    auto ts_xyyy_y = pbuffer.data(idx_ovl_gp + 19);

    auto ts_xzzz_z = pbuffer.data(idx_ovl_gp + 29);

    auto ts_yyyy_x = pbuffer.data(idx_ovl_gp + 30);

    auto ts_yyyy_y = pbuffer.data(idx_ovl_gp + 31);

    auto ts_yyyy_z = pbuffer.data(idx_ovl_gp + 32);

    auto ts_yyyz_z = pbuffer.data(idx_ovl_gp + 35);

    auto ts_yyzz_x = pbuffer.data(idx_ovl_gp + 36);

    auto ts_yyzz_y = pbuffer.data(idx_ovl_gp + 37);

    auto ts_yyzz_z = pbuffer.data(idx_ovl_gp + 38);

    auto ts_yzzz_y = pbuffer.data(idx_ovl_gp + 40);

    auto ts_yzzz_z = pbuffer.data(idx_ovl_gp + 41);

    auto ts_zzzz_x = pbuffer.data(idx_ovl_gp + 42);

    auto ts_zzzz_y = pbuffer.data(idx_ovl_gp + 43);

    auto ts_zzzz_z = pbuffer.data(idx_ovl_gp + 44);

    // Set up components of auxiliary buffer : GD

    auto ts_xxxx_xx = pbuffer.data(idx_ovl_gd);

    auto ts_xxxx_xy = pbuffer.data(idx_ovl_gd + 1);

    auto ts_xxxx_xz = pbuffer.data(idx_ovl_gd + 2);

    auto ts_xxxx_yy = pbuffer.data(idx_ovl_gd + 3);

    auto ts_xxxx_yz = pbuffer.data(idx_ovl_gd + 4);

    auto ts_xxxx_zz = pbuffer.data(idx_ovl_gd + 5);

    auto ts_xxxy_xx = pbuffer.data(idx_ovl_gd + 6);

    auto ts_xxxy_xy = pbuffer.data(idx_ovl_gd + 7);

    auto ts_xxxy_xz = pbuffer.data(idx_ovl_gd + 8);

    auto ts_xxxy_yy = pbuffer.data(idx_ovl_gd + 9);

    auto ts_xxxy_yz = pbuffer.data(idx_ovl_gd + 10);

    auto ts_xxxz_xx = pbuffer.data(idx_ovl_gd + 12);

    auto ts_xxxz_xy = pbuffer.data(idx_ovl_gd + 13);

    auto ts_xxxz_xz = pbuffer.data(idx_ovl_gd + 14);

    auto ts_xxxz_yz = pbuffer.data(idx_ovl_gd + 16);

    auto ts_xxxz_zz = pbuffer.data(idx_ovl_gd + 17);

    auto ts_xxyy_xx = pbuffer.data(idx_ovl_gd + 18);

    auto ts_xxyy_xy = pbuffer.data(idx_ovl_gd + 19);

    auto ts_xxyy_xz = pbuffer.data(idx_ovl_gd + 20);

    auto ts_xxyy_yy = pbuffer.data(idx_ovl_gd + 21);

    auto ts_xxyy_yz = pbuffer.data(idx_ovl_gd + 22);

    auto ts_xxyy_zz = pbuffer.data(idx_ovl_gd + 23);

    auto ts_xxyz_xz = pbuffer.data(idx_ovl_gd + 26);

    auto ts_xxyz_yz = pbuffer.data(idx_ovl_gd + 28);

    auto ts_xxzz_xx = pbuffer.data(idx_ovl_gd + 30);

    auto ts_xxzz_xy = pbuffer.data(idx_ovl_gd + 31);

    auto ts_xxzz_xz = pbuffer.data(idx_ovl_gd + 32);

    auto ts_xxzz_yy = pbuffer.data(idx_ovl_gd + 33);

    auto ts_xxzz_yz = pbuffer.data(idx_ovl_gd + 34);

    auto ts_xxzz_zz = pbuffer.data(idx_ovl_gd + 35);

    auto ts_xyyy_xx = pbuffer.data(idx_ovl_gd + 36);

    auto ts_xyyy_xy = pbuffer.data(idx_ovl_gd + 37);

    auto ts_xyyy_yy = pbuffer.data(idx_ovl_gd + 39);

    auto ts_xyyy_yz = pbuffer.data(idx_ovl_gd + 40);

    auto ts_xyyy_zz = pbuffer.data(idx_ovl_gd + 41);

    auto ts_xyyz_yz = pbuffer.data(idx_ovl_gd + 46);

    auto ts_xyyz_zz = pbuffer.data(idx_ovl_gd + 47);

    auto ts_xyzz_yy = pbuffer.data(idx_ovl_gd + 51);

    auto ts_xyzz_yz = pbuffer.data(idx_ovl_gd + 52);

    auto ts_xzzz_xx = pbuffer.data(idx_ovl_gd + 54);

    auto ts_xzzz_xz = pbuffer.data(idx_ovl_gd + 56);

    auto ts_xzzz_yy = pbuffer.data(idx_ovl_gd + 57);

    auto ts_xzzz_yz = pbuffer.data(idx_ovl_gd + 58);

    auto ts_xzzz_zz = pbuffer.data(idx_ovl_gd + 59);

    auto ts_yyyy_xx = pbuffer.data(idx_ovl_gd + 60);

    auto ts_yyyy_xy = pbuffer.data(idx_ovl_gd + 61);

    auto ts_yyyy_xz = pbuffer.data(idx_ovl_gd + 62);

    auto ts_yyyy_yy = pbuffer.data(idx_ovl_gd + 63);

    auto ts_yyyy_yz = pbuffer.data(idx_ovl_gd + 64);

    auto ts_yyyy_zz = pbuffer.data(idx_ovl_gd + 65);

    auto ts_yyyz_xy = pbuffer.data(idx_ovl_gd + 67);

    auto ts_yyyz_xz = pbuffer.data(idx_ovl_gd + 68);

    auto ts_yyyz_yy = pbuffer.data(idx_ovl_gd + 69);

    auto ts_yyyz_yz = pbuffer.data(idx_ovl_gd + 70);

    auto ts_yyyz_zz = pbuffer.data(idx_ovl_gd + 71);

    auto ts_yyzz_xx = pbuffer.data(idx_ovl_gd + 72);

    auto ts_yyzz_xy = pbuffer.data(idx_ovl_gd + 73);

    auto ts_yyzz_xz = pbuffer.data(idx_ovl_gd + 74);

    auto ts_yyzz_yy = pbuffer.data(idx_ovl_gd + 75);

    auto ts_yyzz_yz = pbuffer.data(idx_ovl_gd + 76);

    auto ts_yyzz_zz = pbuffer.data(idx_ovl_gd + 77);

    auto ts_yzzz_xx = pbuffer.data(idx_ovl_gd + 78);

    auto ts_yzzz_xy = pbuffer.data(idx_ovl_gd + 79);

    auto ts_yzzz_xz = pbuffer.data(idx_ovl_gd + 80);

    auto ts_yzzz_yy = pbuffer.data(idx_ovl_gd + 81);

    auto ts_yzzz_yz = pbuffer.data(idx_ovl_gd + 82);

    auto ts_yzzz_zz = pbuffer.data(idx_ovl_gd + 83);

    auto ts_zzzz_xx = pbuffer.data(idx_ovl_gd + 84);

    auto ts_zzzz_xy = pbuffer.data(idx_ovl_gd + 85);

    auto ts_zzzz_xz = pbuffer.data(idx_ovl_gd + 86);

    auto ts_zzzz_yy = pbuffer.data(idx_ovl_gd + 87);

    auto ts_zzzz_yz = pbuffer.data(idx_ovl_gd + 88);

    auto ts_zzzz_zz = pbuffer.data(idx_ovl_gd + 89);

    // Set up 0-6 components of targeted buffer : HD

    auto ts_xxxxx_xx = pbuffer.data(idx_ovl_hd);

    auto ts_xxxxx_xy = pbuffer.data(idx_ovl_hd + 1);

    auto ts_xxxxx_xz = pbuffer.data(idx_ovl_hd + 2);

    auto ts_xxxxx_yy = pbuffer.data(idx_ovl_hd + 3);

    auto ts_xxxxx_yz = pbuffer.data(idx_ovl_hd + 4);

    auto ts_xxxxx_zz = pbuffer.data(idx_ovl_hd + 5);

    #pragma omp simd aligned(pa_x, ts_xxx_xx, ts_xxx_xy, ts_xxx_xz, ts_xxx_yy, ts_xxx_yz, ts_xxx_zz, ts_xxxx_x, ts_xxxx_xx, ts_xxxx_xy, ts_xxxx_xz, ts_xxxx_y, ts_xxxx_yy, ts_xxxx_yz, ts_xxxx_z, ts_xxxx_zz, ts_xxxxx_xx, ts_xxxxx_xy, ts_xxxxx_xz, ts_xxxxx_yy, ts_xxxxx_yz, ts_xxxxx_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxx_xx[i] = 4.0 * ts_xxx_xx[i] * fe_0 + 2.0 * ts_xxxx_x[i] * fe_0 + ts_xxxx_xx[i] * pa_x[i];

        ts_xxxxx_xy[i] = 4.0 * ts_xxx_xy[i] * fe_0 + ts_xxxx_y[i] * fe_0 + ts_xxxx_xy[i] * pa_x[i];

        ts_xxxxx_xz[i] = 4.0 * ts_xxx_xz[i] * fe_0 + ts_xxxx_z[i] * fe_0 + ts_xxxx_xz[i] * pa_x[i];

        ts_xxxxx_yy[i] = 4.0 * ts_xxx_yy[i] * fe_0 + ts_xxxx_yy[i] * pa_x[i];

        ts_xxxxx_yz[i] = 4.0 * ts_xxx_yz[i] * fe_0 + ts_xxxx_yz[i] * pa_x[i];

        ts_xxxxx_zz[i] = 4.0 * ts_xxx_zz[i] * fe_0 + ts_xxxx_zz[i] * pa_x[i];
    }

    // Set up 6-12 components of targeted buffer : HD

    auto ts_xxxxy_xx = pbuffer.data(idx_ovl_hd + 6);

    auto ts_xxxxy_xy = pbuffer.data(idx_ovl_hd + 7);

    auto ts_xxxxy_xz = pbuffer.data(idx_ovl_hd + 8);

    auto ts_xxxxy_yy = pbuffer.data(idx_ovl_hd + 9);

    auto ts_xxxxy_yz = pbuffer.data(idx_ovl_hd + 10);

    auto ts_xxxxy_zz = pbuffer.data(idx_ovl_hd + 11);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxxx_x, ts_xxxx_xx, ts_xxxx_xy, ts_xxxx_xz, ts_xxxx_zz, ts_xxxxy_xx, ts_xxxxy_xy, ts_xxxxy_xz, ts_xxxxy_yy, ts_xxxxy_yz, ts_xxxxy_zz, ts_xxxy_yy, ts_xxxy_yz, ts_xxy_yy, ts_xxy_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxy_xx[i] = ts_xxxx_xx[i] * pa_y[i];

        ts_xxxxy_xy[i] = ts_xxxx_x[i] * fe_0 + ts_xxxx_xy[i] * pa_y[i];

        ts_xxxxy_xz[i] = ts_xxxx_xz[i] * pa_y[i];

        ts_xxxxy_yy[i] = 3.0 * ts_xxy_yy[i] * fe_0 + ts_xxxy_yy[i] * pa_x[i];

        ts_xxxxy_yz[i] = 3.0 * ts_xxy_yz[i] * fe_0 + ts_xxxy_yz[i] * pa_x[i];

        ts_xxxxy_zz[i] = ts_xxxx_zz[i] * pa_y[i];
    }

    // Set up 12-18 components of targeted buffer : HD

    auto ts_xxxxz_xx = pbuffer.data(idx_ovl_hd + 12);

    auto ts_xxxxz_xy = pbuffer.data(idx_ovl_hd + 13);

    auto ts_xxxxz_xz = pbuffer.data(idx_ovl_hd + 14);

    auto ts_xxxxz_yy = pbuffer.data(idx_ovl_hd + 15);

    auto ts_xxxxz_yz = pbuffer.data(idx_ovl_hd + 16);

    auto ts_xxxxz_zz = pbuffer.data(idx_ovl_hd + 17);

    #pragma omp simd aligned(pa_x, pa_z, ts_xxxx_x, ts_xxxx_xx, ts_xxxx_xy, ts_xxxx_xz, ts_xxxx_yy, ts_xxxxz_xx, ts_xxxxz_xy, ts_xxxxz_xz, ts_xxxxz_yy, ts_xxxxz_yz, ts_xxxxz_zz, ts_xxxz_yz, ts_xxxz_zz, ts_xxz_yz, ts_xxz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxz_xx[i] = ts_xxxx_xx[i] * pa_z[i];

        ts_xxxxz_xy[i] = ts_xxxx_xy[i] * pa_z[i];

        ts_xxxxz_xz[i] = ts_xxxx_x[i] * fe_0 + ts_xxxx_xz[i] * pa_z[i];

        ts_xxxxz_yy[i] = ts_xxxx_yy[i] * pa_z[i];

        ts_xxxxz_yz[i] = 3.0 * ts_xxz_yz[i] * fe_0 + ts_xxxz_yz[i] * pa_x[i];

        ts_xxxxz_zz[i] = 3.0 * ts_xxz_zz[i] * fe_0 + ts_xxxz_zz[i] * pa_x[i];
    }

    // Set up 18-24 components of targeted buffer : HD

    auto ts_xxxyy_xx = pbuffer.data(idx_ovl_hd + 18);

    auto ts_xxxyy_xy = pbuffer.data(idx_ovl_hd + 19);

    auto ts_xxxyy_xz = pbuffer.data(idx_ovl_hd + 20);

    auto ts_xxxyy_yy = pbuffer.data(idx_ovl_hd + 21);

    auto ts_xxxyy_yz = pbuffer.data(idx_ovl_hd + 22);

    auto ts_xxxyy_zz = pbuffer.data(idx_ovl_hd + 23);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxx_xx, ts_xxx_xz, ts_xxxy_xx, ts_xxxy_xz, ts_xxxyy_xx, ts_xxxyy_xy, ts_xxxyy_xz, ts_xxxyy_yy, ts_xxxyy_yz, ts_xxxyy_zz, ts_xxyy_xy, ts_xxyy_y, ts_xxyy_yy, ts_xxyy_yz, ts_xxyy_zz, ts_xyy_xy, ts_xyy_yy, ts_xyy_yz, ts_xyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxyy_xx[i] = ts_xxx_xx[i] * fe_0 + ts_xxxy_xx[i] * pa_y[i];

        ts_xxxyy_xy[i] = 2.0 * ts_xyy_xy[i] * fe_0 + ts_xxyy_y[i] * fe_0 + ts_xxyy_xy[i] * pa_x[i];

        ts_xxxyy_xz[i] = ts_xxx_xz[i] * fe_0 + ts_xxxy_xz[i] * pa_y[i];

        ts_xxxyy_yy[i] = 2.0 * ts_xyy_yy[i] * fe_0 + ts_xxyy_yy[i] * pa_x[i];

        ts_xxxyy_yz[i] = 2.0 * ts_xyy_yz[i] * fe_0 + ts_xxyy_yz[i] * pa_x[i];

        ts_xxxyy_zz[i] = 2.0 * ts_xyy_zz[i] * fe_0 + ts_xxyy_zz[i] * pa_x[i];
    }

    // Set up 24-30 components of targeted buffer : HD

    auto ts_xxxyz_xx = pbuffer.data(idx_ovl_hd + 24);

    auto ts_xxxyz_xy = pbuffer.data(idx_ovl_hd + 25);

    auto ts_xxxyz_xz = pbuffer.data(idx_ovl_hd + 26);

    auto ts_xxxyz_yy = pbuffer.data(idx_ovl_hd + 27);

    auto ts_xxxyz_yz = pbuffer.data(idx_ovl_hd + 28);

    auto ts_xxxyz_zz = pbuffer.data(idx_ovl_hd + 29);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, ts_xxxy_xy, ts_xxxy_yy, ts_xxxyz_xx, ts_xxxyz_xy, ts_xxxyz_xz, ts_xxxyz_yy, ts_xxxyz_yz, ts_xxxyz_zz, ts_xxxz_xx, ts_xxxz_xz, ts_xxxz_zz, ts_xxyz_yz, ts_xyz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxyz_xx[i] = ts_xxxz_xx[i] * pa_y[i];

        ts_xxxyz_xy[i] = ts_xxxy_xy[i] * pa_z[i];

        ts_xxxyz_xz[i] = ts_xxxz_xz[i] * pa_y[i];

        ts_xxxyz_yy[i] = ts_xxxy_yy[i] * pa_z[i];

        ts_xxxyz_yz[i] = 2.0 * ts_xyz_yz[i] * fe_0 + ts_xxyz_yz[i] * pa_x[i];

        ts_xxxyz_zz[i] = ts_xxxz_zz[i] * pa_y[i];
    }

    // Set up 30-36 components of targeted buffer : HD

    auto ts_xxxzz_xx = pbuffer.data(idx_ovl_hd + 30);

    auto ts_xxxzz_xy = pbuffer.data(idx_ovl_hd + 31);

    auto ts_xxxzz_xz = pbuffer.data(idx_ovl_hd + 32);

    auto ts_xxxzz_yy = pbuffer.data(idx_ovl_hd + 33);

    auto ts_xxxzz_yz = pbuffer.data(idx_ovl_hd + 34);

    auto ts_xxxzz_zz = pbuffer.data(idx_ovl_hd + 35);

    #pragma omp simd aligned(pa_x, pa_z, ts_xxx_xx, ts_xxx_xy, ts_xxxz_xx, ts_xxxz_xy, ts_xxxzz_xx, ts_xxxzz_xy, ts_xxxzz_xz, ts_xxxzz_yy, ts_xxxzz_yz, ts_xxxzz_zz, ts_xxzz_xz, ts_xxzz_yy, ts_xxzz_yz, ts_xxzz_z, ts_xxzz_zz, ts_xzz_xz, ts_xzz_yy, ts_xzz_yz, ts_xzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxzz_xx[i] = ts_xxx_xx[i] * fe_0 + ts_xxxz_xx[i] * pa_z[i];

        ts_xxxzz_xy[i] = ts_xxx_xy[i] * fe_0 + ts_xxxz_xy[i] * pa_z[i];

        ts_xxxzz_xz[i] = 2.0 * ts_xzz_xz[i] * fe_0 + ts_xxzz_z[i] * fe_0 + ts_xxzz_xz[i] * pa_x[i];

        ts_xxxzz_yy[i] = 2.0 * ts_xzz_yy[i] * fe_0 + ts_xxzz_yy[i] * pa_x[i];

        ts_xxxzz_yz[i] = 2.0 * ts_xzz_yz[i] * fe_0 + ts_xxzz_yz[i] * pa_x[i];

        ts_xxxzz_zz[i] = 2.0 * ts_xzz_zz[i] * fe_0 + ts_xxzz_zz[i] * pa_x[i];
    }

    // Set up 36-42 components of targeted buffer : HD

    auto ts_xxyyy_xx = pbuffer.data(idx_ovl_hd + 36);

    auto ts_xxyyy_xy = pbuffer.data(idx_ovl_hd + 37);

    auto ts_xxyyy_xz = pbuffer.data(idx_ovl_hd + 38);

    auto ts_xxyyy_yy = pbuffer.data(idx_ovl_hd + 39);

    auto ts_xxyyy_yz = pbuffer.data(idx_ovl_hd + 40);

    auto ts_xxyyy_zz = pbuffer.data(idx_ovl_hd + 41);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxy_xx, ts_xxy_xz, ts_xxyy_xx, ts_xxyy_xz, ts_xxyyy_xx, ts_xxyyy_xy, ts_xxyyy_xz, ts_xxyyy_yy, ts_xxyyy_yz, ts_xxyyy_zz, ts_xyyy_xy, ts_xyyy_y, ts_xyyy_yy, ts_xyyy_yz, ts_xyyy_zz, ts_yyy_xy, ts_yyy_yy, ts_yyy_yz, ts_yyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyyy_xx[i] = 2.0 * ts_xxy_xx[i] * fe_0 + ts_xxyy_xx[i] * pa_y[i];

        ts_xxyyy_xy[i] = ts_yyy_xy[i] * fe_0 + ts_xyyy_y[i] * fe_0 + ts_xyyy_xy[i] * pa_x[i];

        ts_xxyyy_xz[i] = 2.0 * ts_xxy_xz[i] * fe_0 + ts_xxyy_xz[i] * pa_y[i];

        ts_xxyyy_yy[i] = ts_yyy_yy[i] * fe_0 + ts_xyyy_yy[i] * pa_x[i];

        ts_xxyyy_yz[i] = ts_yyy_yz[i] * fe_0 + ts_xyyy_yz[i] * pa_x[i];

        ts_xxyyy_zz[i] = ts_yyy_zz[i] * fe_0 + ts_xyyy_zz[i] * pa_x[i];
    }

    // Set up 42-48 components of targeted buffer : HD

    auto ts_xxyyz_xx = pbuffer.data(idx_ovl_hd + 42);

    auto ts_xxyyz_xy = pbuffer.data(idx_ovl_hd + 43);

    auto ts_xxyyz_xz = pbuffer.data(idx_ovl_hd + 44);

    auto ts_xxyyz_yy = pbuffer.data(idx_ovl_hd + 45);

    auto ts_xxyyz_yz = pbuffer.data(idx_ovl_hd + 46);

    auto ts_xxyyz_zz = pbuffer.data(idx_ovl_hd + 47);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, ts_xxyy_xx, ts_xxyy_xy, ts_xxyy_yy, ts_xxyyz_xx, ts_xxyyz_xy, ts_xxyyz_xz, ts_xxyyz_yy, ts_xxyyz_yz, ts_xxyyz_zz, ts_xxyz_xz, ts_xxz_xz, ts_xyyz_yz, ts_xyyz_zz, ts_yyz_yz, ts_yyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyyz_xx[i] = ts_xxyy_xx[i] * pa_z[i];

        ts_xxyyz_xy[i] = ts_xxyy_xy[i] * pa_z[i];

        ts_xxyyz_xz[i] = ts_xxz_xz[i] * fe_0 + ts_xxyz_xz[i] * pa_y[i];

        ts_xxyyz_yy[i] = ts_xxyy_yy[i] * pa_z[i];

        ts_xxyyz_yz[i] = ts_yyz_yz[i] * fe_0 + ts_xyyz_yz[i] * pa_x[i];

        ts_xxyyz_zz[i] = ts_yyz_zz[i] * fe_0 + ts_xyyz_zz[i] * pa_x[i];
    }

    // Set up 48-54 components of targeted buffer : HD

    auto ts_xxyzz_xx = pbuffer.data(idx_ovl_hd + 48);

    auto ts_xxyzz_xy = pbuffer.data(idx_ovl_hd + 49);

    auto ts_xxyzz_xz = pbuffer.data(idx_ovl_hd + 50);

    auto ts_xxyzz_yy = pbuffer.data(idx_ovl_hd + 51);

    auto ts_xxyzz_yz = pbuffer.data(idx_ovl_hd + 52);

    auto ts_xxyzz_zz = pbuffer.data(idx_ovl_hd + 53);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxyzz_xx, ts_xxyzz_xy, ts_xxyzz_xz, ts_xxyzz_yy, ts_xxyzz_yz, ts_xxyzz_zz, ts_xxzz_x, ts_xxzz_xx, ts_xxzz_xy, ts_xxzz_xz, ts_xxzz_zz, ts_xyzz_yy, ts_xyzz_yz, ts_yzz_yy, ts_yzz_yz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyzz_xx[i] = ts_xxzz_xx[i] * pa_y[i];

        ts_xxyzz_xy[i] = ts_xxzz_x[i] * fe_0 + ts_xxzz_xy[i] * pa_y[i];

        ts_xxyzz_xz[i] = ts_xxzz_xz[i] * pa_y[i];

        ts_xxyzz_yy[i] = ts_yzz_yy[i] * fe_0 + ts_xyzz_yy[i] * pa_x[i];

        ts_xxyzz_yz[i] = ts_yzz_yz[i] * fe_0 + ts_xyzz_yz[i] * pa_x[i];

        ts_xxyzz_zz[i] = ts_xxzz_zz[i] * pa_y[i];
    }

    // Set up 54-60 components of targeted buffer : HD

    auto ts_xxzzz_xx = pbuffer.data(idx_ovl_hd + 54);

    auto ts_xxzzz_xy = pbuffer.data(idx_ovl_hd + 55);

    auto ts_xxzzz_xz = pbuffer.data(idx_ovl_hd + 56);

    auto ts_xxzzz_yy = pbuffer.data(idx_ovl_hd + 57);

    auto ts_xxzzz_yz = pbuffer.data(idx_ovl_hd + 58);

    auto ts_xxzzz_zz = pbuffer.data(idx_ovl_hd + 59);

    #pragma omp simd aligned(pa_x, pa_z, ts_xxz_xx, ts_xxz_xy, ts_xxzz_xx, ts_xxzz_xy, ts_xxzzz_xx, ts_xxzzz_xy, ts_xxzzz_xz, ts_xxzzz_yy, ts_xxzzz_yz, ts_xxzzz_zz, ts_xzzz_xz, ts_xzzz_yy, ts_xzzz_yz, ts_xzzz_z, ts_xzzz_zz, ts_zzz_xz, ts_zzz_yy, ts_zzz_yz, ts_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxzzz_xx[i] = 2.0 * ts_xxz_xx[i] * fe_0 + ts_xxzz_xx[i] * pa_z[i];

        ts_xxzzz_xy[i] = 2.0 * ts_xxz_xy[i] * fe_0 + ts_xxzz_xy[i] * pa_z[i];

        ts_xxzzz_xz[i] = ts_zzz_xz[i] * fe_0 + ts_xzzz_z[i] * fe_0 + ts_xzzz_xz[i] * pa_x[i];

        ts_xxzzz_yy[i] = ts_zzz_yy[i] * fe_0 + ts_xzzz_yy[i] * pa_x[i];

        ts_xxzzz_yz[i] = ts_zzz_yz[i] * fe_0 + ts_xzzz_yz[i] * pa_x[i];

        ts_xxzzz_zz[i] = ts_zzz_zz[i] * fe_0 + ts_xzzz_zz[i] * pa_x[i];
    }

    // Set up 60-66 components of targeted buffer : HD

    auto ts_xyyyy_xx = pbuffer.data(idx_ovl_hd + 60);

    auto ts_xyyyy_xy = pbuffer.data(idx_ovl_hd + 61);

    auto ts_xyyyy_xz = pbuffer.data(idx_ovl_hd + 62);

    auto ts_xyyyy_yy = pbuffer.data(idx_ovl_hd + 63);

    auto ts_xyyyy_yz = pbuffer.data(idx_ovl_hd + 64);

    auto ts_xyyyy_zz = pbuffer.data(idx_ovl_hd + 65);

    #pragma omp simd aligned(pa_x, ts_xyyyy_xx, ts_xyyyy_xy, ts_xyyyy_xz, ts_xyyyy_yy, ts_xyyyy_yz, ts_xyyyy_zz, ts_yyyy_x, ts_yyyy_xx, ts_yyyy_xy, ts_yyyy_xz, ts_yyyy_y, ts_yyyy_yy, ts_yyyy_yz, ts_yyyy_z, ts_yyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyyy_xx[i] = 2.0 * ts_yyyy_x[i] * fe_0 + ts_yyyy_xx[i] * pa_x[i];

        ts_xyyyy_xy[i] = ts_yyyy_y[i] * fe_0 + ts_yyyy_xy[i] * pa_x[i];

        ts_xyyyy_xz[i] = ts_yyyy_z[i] * fe_0 + ts_yyyy_xz[i] * pa_x[i];

        ts_xyyyy_yy[i] = ts_yyyy_yy[i] * pa_x[i];

        ts_xyyyy_yz[i] = ts_yyyy_yz[i] * pa_x[i];

        ts_xyyyy_zz[i] = ts_yyyy_zz[i] * pa_x[i];
    }

    // Set up 66-72 components of targeted buffer : HD

    auto ts_xyyyz_xx = pbuffer.data(idx_ovl_hd + 66);

    auto ts_xyyyz_xy = pbuffer.data(idx_ovl_hd + 67);

    auto ts_xyyyz_xz = pbuffer.data(idx_ovl_hd + 68);

    auto ts_xyyyz_yy = pbuffer.data(idx_ovl_hd + 69);

    auto ts_xyyyz_yz = pbuffer.data(idx_ovl_hd + 70);

    auto ts_xyyyz_zz = pbuffer.data(idx_ovl_hd + 71);

    #pragma omp simd aligned(pa_x, pa_z, ts_xyyy_xx, ts_xyyy_xy, ts_xyyyz_xx, ts_xyyyz_xy, ts_xyyyz_xz, ts_xyyyz_yy, ts_xyyyz_yz, ts_xyyyz_zz, ts_yyyz_xz, ts_yyyz_yy, ts_yyyz_yz, ts_yyyz_z, ts_yyyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyyz_xx[i] = ts_xyyy_xx[i] * pa_z[i];

        ts_xyyyz_xy[i] = ts_xyyy_xy[i] * pa_z[i];

        ts_xyyyz_xz[i] = ts_yyyz_z[i] * fe_0 + ts_yyyz_xz[i] * pa_x[i];

        ts_xyyyz_yy[i] = ts_yyyz_yy[i] * pa_x[i];

        ts_xyyyz_yz[i] = ts_yyyz_yz[i] * pa_x[i];

        ts_xyyyz_zz[i] = ts_yyyz_zz[i] * pa_x[i];
    }

    // Set up 72-78 components of targeted buffer : HD

    auto ts_xyyzz_xx = pbuffer.data(idx_ovl_hd + 72);

    auto ts_xyyzz_xy = pbuffer.data(idx_ovl_hd + 73);

    auto ts_xyyzz_xz = pbuffer.data(idx_ovl_hd + 74);

    auto ts_xyyzz_yy = pbuffer.data(idx_ovl_hd + 75);

    auto ts_xyyzz_yz = pbuffer.data(idx_ovl_hd + 76);

    auto ts_xyyzz_zz = pbuffer.data(idx_ovl_hd + 77);

    #pragma omp simd aligned(pa_x, ts_xyyzz_xx, ts_xyyzz_xy, ts_xyyzz_xz, ts_xyyzz_yy, ts_xyyzz_yz, ts_xyyzz_zz, ts_yyzz_x, ts_yyzz_xx, ts_yyzz_xy, ts_yyzz_xz, ts_yyzz_y, ts_yyzz_yy, ts_yyzz_yz, ts_yyzz_z, ts_yyzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyzz_xx[i] = 2.0 * ts_yyzz_x[i] * fe_0 + ts_yyzz_xx[i] * pa_x[i];

        ts_xyyzz_xy[i] = ts_yyzz_y[i] * fe_0 + ts_yyzz_xy[i] * pa_x[i];

        ts_xyyzz_xz[i] = ts_yyzz_z[i] * fe_0 + ts_yyzz_xz[i] * pa_x[i];

        ts_xyyzz_yy[i] = ts_yyzz_yy[i] * pa_x[i];

        ts_xyyzz_yz[i] = ts_yyzz_yz[i] * pa_x[i];

        ts_xyyzz_zz[i] = ts_yyzz_zz[i] * pa_x[i];
    }

    // Set up 78-84 components of targeted buffer : HD

    auto ts_xyzzz_xx = pbuffer.data(idx_ovl_hd + 78);

    auto ts_xyzzz_xy = pbuffer.data(idx_ovl_hd + 79);

    auto ts_xyzzz_xz = pbuffer.data(idx_ovl_hd + 80);

    auto ts_xyzzz_yy = pbuffer.data(idx_ovl_hd + 81);

    auto ts_xyzzz_yz = pbuffer.data(idx_ovl_hd + 82);

    auto ts_xyzzz_zz = pbuffer.data(idx_ovl_hd + 83);

    #pragma omp simd aligned(pa_x, pa_y, ts_xyzzz_xx, ts_xyzzz_xy, ts_xyzzz_xz, ts_xyzzz_yy, ts_xyzzz_yz, ts_xyzzz_zz, ts_xzzz_xx, ts_xzzz_xz, ts_yzzz_xy, ts_yzzz_y, ts_yzzz_yy, ts_yzzz_yz, ts_yzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyzzz_xx[i] = ts_xzzz_xx[i] * pa_y[i];

        ts_xyzzz_xy[i] = ts_yzzz_y[i] * fe_0 + ts_yzzz_xy[i] * pa_x[i];

        ts_xyzzz_xz[i] = ts_xzzz_xz[i] * pa_y[i];

        ts_xyzzz_yy[i] = ts_yzzz_yy[i] * pa_x[i];

        ts_xyzzz_yz[i] = ts_yzzz_yz[i] * pa_x[i];

        ts_xyzzz_zz[i] = ts_yzzz_zz[i] * pa_x[i];
    }

    // Set up 84-90 components of targeted buffer : HD

    auto ts_xzzzz_xx = pbuffer.data(idx_ovl_hd + 84);

    auto ts_xzzzz_xy = pbuffer.data(idx_ovl_hd + 85);

    auto ts_xzzzz_xz = pbuffer.data(idx_ovl_hd + 86);

    auto ts_xzzzz_yy = pbuffer.data(idx_ovl_hd + 87);

    auto ts_xzzzz_yz = pbuffer.data(idx_ovl_hd + 88);

    auto ts_xzzzz_zz = pbuffer.data(idx_ovl_hd + 89);

    #pragma omp simd aligned(pa_x, ts_xzzzz_xx, ts_xzzzz_xy, ts_xzzzz_xz, ts_xzzzz_yy, ts_xzzzz_yz, ts_xzzzz_zz, ts_zzzz_x, ts_zzzz_xx, ts_zzzz_xy, ts_zzzz_xz, ts_zzzz_y, ts_zzzz_yy, ts_zzzz_yz, ts_zzzz_z, ts_zzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xzzzz_xx[i] = 2.0 * ts_zzzz_x[i] * fe_0 + ts_zzzz_xx[i] * pa_x[i];

        ts_xzzzz_xy[i] = ts_zzzz_y[i] * fe_0 + ts_zzzz_xy[i] * pa_x[i];

        ts_xzzzz_xz[i] = ts_zzzz_z[i] * fe_0 + ts_zzzz_xz[i] * pa_x[i];

        ts_xzzzz_yy[i] = ts_zzzz_yy[i] * pa_x[i];

        ts_xzzzz_yz[i] = ts_zzzz_yz[i] * pa_x[i];

        ts_xzzzz_zz[i] = ts_zzzz_zz[i] * pa_x[i];
    }

    // Set up 90-96 components of targeted buffer : HD

    auto ts_yyyyy_xx = pbuffer.data(idx_ovl_hd + 90);

    auto ts_yyyyy_xy = pbuffer.data(idx_ovl_hd + 91);

    auto ts_yyyyy_xz = pbuffer.data(idx_ovl_hd + 92);

    auto ts_yyyyy_yy = pbuffer.data(idx_ovl_hd + 93);

    auto ts_yyyyy_yz = pbuffer.data(idx_ovl_hd + 94);

    auto ts_yyyyy_zz = pbuffer.data(idx_ovl_hd + 95);

    #pragma omp simd aligned(pa_y, ts_yyy_xx, ts_yyy_xy, ts_yyy_xz, ts_yyy_yy, ts_yyy_yz, ts_yyy_zz, ts_yyyy_x, ts_yyyy_xx, ts_yyyy_xy, ts_yyyy_xz, ts_yyyy_y, ts_yyyy_yy, ts_yyyy_yz, ts_yyyy_z, ts_yyyy_zz, ts_yyyyy_xx, ts_yyyyy_xy, ts_yyyyy_xz, ts_yyyyy_yy, ts_yyyyy_yz, ts_yyyyy_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyyy_xx[i] = 4.0 * ts_yyy_xx[i] * fe_0 + ts_yyyy_xx[i] * pa_y[i];

        ts_yyyyy_xy[i] = 4.0 * ts_yyy_xy[i] * fe_0 + ts_yyyy_x[i] * fe_0 + ts_yyyy_xy[i] * pa_y[i];

        ts_yyyyy_xz[i] = 4.0 * ts_yyy_xz[i] * fe_0 + ts_yyyy_xz[i] * pa_y[i];

        ts_yyyyy_yy[i] = 4.0 * ts_yyy_yy[i] * fe_0 + 2.0 * ts_yyyy_y[i] * fe_0 + ts_yyyy_yy[i] * pa_y[i];

        ts_yyyyy_yz[i] = 4.0 * ts_yyy_yz[i] * fe_0 + ts_yyyy_z[i] * fe_0 + ts_yyyy_yz[i] * pa_y[i];

        ts_yyyyy_zz[i] = 4.0 * ts_yyy_zz[i] * fe_0 + ts_yyyy_zz[i] * pa_y[i];
    }

    // Set up 96-102 components of targeted buffer : HD

    auto ts_yyyyz_xx = pbuffer.data(idx_ovl_hd + 96);

    auto ts_yyyyz_xy = pbuffer.data(idx_ovl_hd + 97);

    auto ts_yyyyz_xz = pbuffer.data(idx_ovl_hd + 98);

    auto ts_yyyyz_yy = pbuffer.data(idx_ovl_hd + 99);

    auto ts_yyyyz_yz = pbuffer.data(idx_ovl_hd + 100);

    auto ts_yyyyz_zz = pbuffer.data(idx_ovl_hd + 101);

    #pragma omp simd aligned(pa_y, pa_z, ts_yyyy_xx, ts_yyyy_xy, ts_yyyy_y, ts_yyyy_yy, ts_yyyy_yz, ts_yyyyz_xx, ts_yyyyz_xy, ts_yyyyz_xz, ts_yyyyz_yy, ts_yyyyz_yz, ts_yyyyz_zz, ts_yyyz_xz, ts_yyyz_zz, ts_yyz_xz, ts_yyz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyyz_xx[i] = ts_yyyy_xx[i] * pa_z[i];

        ts_yyyyz_xy[i] = ts_yyyy_xy[i] * pa_z[i];

        ts_yyyyz_xz[i] = 3.0 * ts_yyz_xz[i] * fe_0 + ts_yyyz_xz[i] * pa_y[i];

        ts_yyyyz_yy[i] = ts_yyyy_yy[i] * pa_z[i];

        ts_yyyyz_yz[i] = ts_yyyy_y[i] * fe_0 + ts_yyyy_yz[i] * pa_z[i];

        ts_yyyyz_zz[i] = 3.0 * ts_yyz_zz[i] * fe_0 + ts_yyyz_zz[i] * pa_y[i];
    }

    // Set up 102-108 components of targeted buffer : HD

    auto ts_yyyzz_xx = pbuffer.data(idx_ovl_hd + 102);

    auto ts_yyyzz_xy = pbuffer.data(idx_ovl_hd + 103);

    auto ts_yyyzz_xz = pbuffer.data(idx_ovl_hd + 104);

    auto ts_yyyzz_yy = pbuffer.data(idx_ovl_hd + 105);

    auto ts_yyyzz_yz = pbuffer.data(idx_ovl_hd + 106);

    auto ts_yyyzz_zz = pbuffer.data(idx_ovl_hd + 107);

    #pragma omp simd aligned(pa_y, pa_z, ts_yyy_xy, ts_yyy_yy, ts_yyyz_xy, ts_yyyz_yy, ts_yyyzz_xx, ts_yyyzz_xy, ts_yyyzz_xz, ts_yyyzz_yy, ts_yyyzz_yz, ts_yyyzz_zz, ts_yyzz_xx, ts_yyzz_xz, ts_yyzz_yz, ts_yyzz_z, ts_yyzz_zz, ts_yzz_xx, ts_yzz_xz, ts_yzz_yz, ts_yzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyzz_xx[i] = 2.0 * ts_yzz_xx[i] * fe_0 + ts_yyzz_xx[i] * pa_y[i];

        ts_yyyzz_xy[i] = ts_yyy_xy[i] * fe_0 + ts_yyyz_xy[i] * pa_z[i];

        ts_yyyzz_xz[i] = 2.0 * ts_yzz_xz[i] * fe_0 + ts_yyzz_xz[i] * pa_y[i];

        ts_yyyzz_yy[i] = ts_yyy_yy[i] * fe_0 + ts_yyyz_yy[i] * pa_z[i];

        ts_yyyzz_yz[i] = 2.0 * ts_yzz_yz[i] * fe_0 + ts_yyzz_z[i] * fe_0 + ts_yyzz_yz[i] * pa_y[i];

        ts_yyyzz_zz[i] = 2.0 * ts_yzz_zz[i] * fe_0 + ts_yyzz_zz[i] * pa_y[i];
    }

    // Set up 108-114 components of targeted buffer : HD

    auto ts_yyzzz_xx = pbuffer.data(idx_ovl_hd + 108);

    auto ts_yyzzz_xy = pbuffer.data(idx_ovl_hd + 109);

    auto ts_yyzzz_xz = pbuffer.data(idx_ovl_hd + 110);

    auto ts_yyzzz_yy = pbuffer.data(idx_ovl_hd + 111);

    auto ts_yyzzz_yz = pbuffer.data(idx_ovl_hd + 112);

    auto ts_yyzzz_zz = pbuffer.data(idx_ovl_hd + 113);

    #pragma omp simd aligned(pa_y, pa_z, ts_yyz_xy, ts_yyz_yy, ts_yyzz_xy, ts_yyzz_yy, ts_yyzzz_xx, ts_yyzzz_xy, ts_yyzzz_xz, ts_yyzzz_yy, ts_yyzzz_yz, ts_yyzzz_zz, ts_yzzz_xx, ts_yzzz_xz, ts_yzzz_yz, ts_yzzz_z, ts_yzzz_zz, ts_zzz_xx, ts_zzz_xz, ts_zzz_yz, ts_zzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyzzz_xx[i] = ts_zzz_xx[i] * fe_0 + ts_yzzz_xx[i] * pa_y[i];

        ts_yyzzz_xy[i] = 2.0 * ts_yyz_xy[i] * fe_0 + ts_yyzz_xy[i] * pa_z[i];

        ts_yyzzz_xz[i] = ts_zzz_xz[i] * fe_0 + ts_yzzz_xz[i] * pa_y[i];

        ts_yyzzz_yy[i] = 2.0 * ts_yyz_yy[i] * fe_0 + ts_yyzz_yy[i] * pa_z[i];

        ts_yyzzz_yz[i] = ts_zzz_yz[i] * fe_0 + ts_yzzz_z[i] * fe_0 + ts_yzzz_yz[i] * pa_y[i];

        ts_yyzzz_zz[i] = ts_zzz_zz[i] * fe_0 + ts_yzzz_zz[i] * pa_y[i];
    }

    // Set up 114-120 components of targeted buffer : HD

    auto ts_yzzzz_xx = pbuffer.data(idx_ovl_hd + 114);

    auto ts_yzzzz_xy = pbuffer.data(idx_ovl_hd + 115);

    auto ts_yzzzz_xz = pbuffer.data(idx_ovl_hd + 116);

    auto ts_yzzzz_yy = pbuffer.data(idx_ovl_hd + 117);

    auto ts_yzzzz_yz = pbuffer.data(idx_ovl_hd + 118);

    auto ts_yzzzz_zz = pbuffer.data(idx_ovl_hd + 119);

    #pragma omp simd aligned(pa_y, ts_yzzzz_xx, ts_yzzzz_xy, ts_yzzzz_xz, ts_yzzzz_yy, ts_yzzzz_yz, ts_yzzzz_zz, ts_zzzz_x, ts_zzzz_xx, ts_zzzz_xy, ts_zzzz_xz, ts_zzzz_y, ts_zzzz_yy, ts_zzzz_yz, ts_zzzz_z, ts_zzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yzzzz_xx[i] = ts_zzzz_xx[i] * pa_y[i];

        ts_yzzzz_xy[i] = ts_zzzz_x[i] * fe_0 + ts_zzzz_xy[i] * pa_y[i];

        ts_yzzzz_xz[i] = ts_zzzz_xz[i] * pa_y[i];

        ts_yzzzz_yy[i] = 2.0 * ts_zzzz_y[i] * fe_0 + ts_zzzz_yy[i] * pa_y[i];

        ts_yzzzz_yz[i] = ts_zzzz_z[i] * fe_0 + ts_zzzz_yz[i] * pa_y[i];

        ts_yzzzz_zz[i] = ts_zzzz_zz[i] * pa_y[i];
    }

    // Set up 120-126 components of targeted buffer : HD

    auto ts_zzzzz_xx = pbuffer.data(idx_ovl_hd + 120);

    auto ts_zzzzz_xy = pbuffer.data(idx_ovl_hd + 121);

    auto ts_zzzzz_xz = pbuffer.data(idx_ovl_hd + 122);

    auto ts_zzzzz_yy = pbuffer.data(idx_ovl_hd + 123);

    auto ts_zzzzz_yz = pbuffer.data(idx_ovl_hd + 124);

    auto ts_zzzzz_zz = pbuffer.data(idx_ovl_hd + 125);

    #pragma omp simd aligned(pa_z, ts_zzz_xx, ts_zzz_xy, ts_zzz_xz, ts_zzz_yy, ts_zzz_yz, ts_zzz_zz, ts_zzzz_x, ts_zzzz_xx, ts_zzzz_xy, ts_zzzz_xz, ts_zzzz_y, ts_zzzz_yy, ts_zzzz_yz, ts_zzzz_z, ts_zzzz_zz, ts_zzzzz_xx, ts_zzzzz_xy, ts_zzzzz_xz, ts_zzzzz_yy, ts_zzzzz_yz, ts_zzzzz_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_zzzzz_xx[i] = 4.0 * ts_zzz_xx[i] * fe_0 + ts_zzzz_xx[i] * pa_z[i];

        ts_zzzzz_xy[i] = 4.0 * ts_zzz_xy[i] * fe_0 + ts_zzzz_xy[i] * pa_z[i];

        ts_zzzzz_xz[i] = 4.0 * ts_zzz_xz[i] * fe_0 + ts_zzzz_x[i] * fe_0 + ts_zzzz_xz[i] * pa_z[i];

        ts_zzzzz_yy[i] = 4.0 * ts_zzz_yy[i] * fe_0 + ts_zzzz_yy[i] * pa_z[i];

        ts_zzzzz_yz[i] = 4.0 * ts_zzz_yz[i] * fe_0 + ts_zzzz_y[i] * fe_0 + ts_zzzz_yz[i] * pa_z[i];

        ts_zzzzz_zz[i] = 4.0 * ts_zzz_zz[i] * fe_0 + 2.0 * ts_zzzz_z[i] * fe_0 + ts_zzzz_zz[i] * pa_z[i];
    }

}

} // ovlrec namespace

