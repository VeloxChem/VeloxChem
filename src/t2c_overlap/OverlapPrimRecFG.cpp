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

#include "OverlapPrimRecFG.hpp"

namespace ovlrec {  // ovlrec namespace

auto
comp_prim_overlap_fg(CSimdArray<double>&       pbuffer,
                     const size_t              idx_ovl_fg,
                     const size_t              idx_ovl_pg,
                     const size_t              idx_ovl_df,
                     const size_t              idx_ovl_dg,
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

    // Set up components of auxiliary buffer : PG

    auto ts_x_xxxx = pbuffer.data(idx_ovl_pg);

    auto ts_x_xxxy = pbuffer.data(idx_ovl_pg + 1);

    auto ts_x_xxxz = pbuffer.data(idx_ovl_pg + 2);

    auto ts_x_xxyy = pbuffer.data(idx_ovl_pg + 3);

    auto ts_x_xxyz = pbuffer.data(idx_ovl_pg + 4);

    auto ts_x_xxzz = pbuffer.data(idx_ovl_pg + 5);

    auto ts_x_xyyy = pbuffer.data(idx_ovl_pg + 6);

    auto ts_x_xyyz = pbuffer.data(idx_ovl_pg + 7);

    auto ts_x_xyzz = pbuffer.data(idx_ovl_pg + 8);

    auto ts_x_xzzz = pbuffer.data(idx_ovl_pg + 9);

    auto ts_x_yyyy = pbuffer.data(idx_ovl_pg + 10);

    auto ts_x_yyyz = pbuffer.data(idx_ovl_pg + 11);

    auto ts_x_yyzz = pbuffer.data(idx_ovl_pg + 12);

    auto ts_x_yzzz = pbuffer.data(idx_ovl_pg + 13);

    auto ts_x_zzzz = pbuffer.data(idx_ovl_pg + 14);

    auto ts_y_xxxx = pbuffer.data(idx_ovl_pg + 15);

    auto ts_y_xxxy = pbuffer.data(idx_ovl_pg + 16);

    auto ts_y_xxxz = pbuffer.data(idx_ovl_pg + 17);

    auto ts_y_xxyy = pbuffer.data(idx_ovl_pg + 18);

    auto ts_y_xxyz = pbuffer.data(idx_ovl_pg + 19);

    auto ts_y_xxzz = pbuffer.data(idx_ovl_pg + 20);

    auto ts_y_xyyy = pbuffer.data(idx_ovl_pg + 21);

    auto ts_y_xyyz = pbuffer.data(idx_ovl_pg + 22);

    auto ts_y_xyzz = pbuffer.data(idx_ovl_pg + 23);

    auto ts_y_xzzz = pbuffer.data(idx_ovl_pg + 24);

    auto ts_y_yyyy = pbuffer.data(idx_ovl_pg + 25);

    auto ts_y_yyyz = pbuffer.data(idx_ovl_pg + 26);

    auto ts_y_yyzz = pbuffer.data(idx_ovl_pg + 27);

    auto ts_y_yzzz = pbuffer.data(idx_ovl_pg + 28);

    auto ts_y_zzzz = pbuffer.data(idx_ovl_pg + 29);

    auto ts_z_xxxx = pbuffer.data(idx_ovl_pg + 30);

    auto ts_z_xxxy = pbuffer.data(idx_ovl_pg + 31);

    auto ts_z_xxxz = pbuffer.data(idx_ovl_pg + 32);

    auto ts_z_xxyy = pbuffer.data(idx_ovl_pg + 33);

    auto ts_z_xxyz = pbuffer.data(idx_ovl_pg + 34);

    auto ts_z_xxzz = pbuffer.data(idx_ovl_pg + 35);

    auto ts_z_xyyy = pbuffer.data(idx_ovl_pg + 36);

    auto ts_z_xyyz = pbuffer.data(idx_ovl_pg + 37);

    auto ts_z_xyzz = pbuffer.data(idx_ovl_pg + 38);

    auto ts_z_xzzz = pbuffer.data(idx_ovl_pg + 39);

    auto ts_z_yyyy = pbuffer.data(idx_ovl_pg + 40);

    auto ts_z_yyyz = pbuffer.data(idx_ovl_pg + 41);

    auto ts_z_yyzz = pbuffer.data(idx_ovl_pg + 42);

    auto ts_z_yzzz = pbuffer.data(idx_ovl_pg + 43);

    auto ts_z_zzzz = pbuffer.data(idx_ovl_pg + 44);

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

    auto ts_yz_xyz = pbuffer.data(idx_ovl_df + 44);

    auto ts_yz_yyz = pbuffer.data(idx_ovl_df + 47);

    auto ts_yz_yzz = pbuffer.data(idx_ovl_df + 48);

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

    // Set up components of auxiliary buffer : DG

    auto ts_xx_xxxx = pbuffer.data(idx_ovl_dg);

    auto ts_xx_xxxy = pbuffer.data(idx_ovl_dg + 1);

    auto ts_xx_xxxz = pbuffer.data(idx_ovl_dg + 2);

    auto ts_xx_xxyy = pbuffer.data(idx_ovl_dg + 3);

    auto ts_xx_xxyz = pbuffer.data(idx_ovl_dg + 4);

    auto ts_xx_xxzz = pbuffer.data(idx_ovl_dg + 5);

    auto ts_xx_xyyy = pbuffer.data(idx_ovl_dg + 6);

    auto ts_xx_xyyz = pbuffer.data(idx_ovl_dg + 7);

    auto ts_xx_xyzz = pbuffer.data(idx_ovl_dg + 8);

    auto ts_xx_xzzz = pbuffer.data(idx_ovl_dg + 9);

    auto ts_xx_yyyy = pbuffer.data(idx_ovl_dg + 10);

    auto ts_xx_yyyz = pbuffer.data(idx_ovl_dg + 11);

    auto ts_xx_yyzz = pbuffer.data(idx_ovl_dg + 12);

    auto ts_xx_yzzz = pbuffer.data(idx_ovl_dg + 13);

    auto ts_xx_zzzz = pbuffer.data(idx_ovl_dg + 14);

    auto ts_xy_xxxy = pbuffer.data(idx_ovl_dg + 16);

    auto ts_xy_xxyy = pbuffer.data(idx_ovl_dg + 18);

    auto ts_xy_xyyy = pbuffer.data(idx_ovl_dg + 21);

    auto ts_xy_yyyy = pbuffer.data(idx_ovl_dg + 25);

    auto ts_xy_yyyz = pbuffer.data(idx_ovl_dg + 26);

    auto ts_xy_yyzz = pbuffer.data(idx_ovl_dg + 27);

    auto ts_xy_yzzz = pbuffer.data(idx_ovl_dg + 28);

    auto ts_xz_xxxx = pbuffer.data(idx_ovl_dg + 30);

    auto ts_xz_xxxz = pbuffer.data(idx_ovl_dg + 32);

    auto ts_xz_xxzz = pbuffer.data(idx_ovl_dg + 35);

    auto ts_xz_xzzz = pbuffer.data(idx_ovl_dg + 39);

    auto ts_xz_yyyz = pbuffer.data(idx_ovl_dg + 41);

    auto ts_xz_yyzz = pbuffer.data(idx_ovl_dg + 42);

    auto ts_xz_yzzz = pbuffer.data(idx_ovl_dg + 43);

    auto ts_xz_zzzz = pbuffer.data(idx_ovl_dg + 44);

    auto ts_yy_xxxx = pbuffer.data(idx_ovl_dg + 45);

    auto ts_yy_xxxy = pbuffer.data(idx_ovl_dg + 46);

    auto ts_yy_xxxz = pbuffer.data(idx_ovl_dg + 47);

    auto ts_yy_xxyy = pbuffer.data(idx_ovl_dg + 48);

    auto ts_yy_xxyz = pbuffer.data(idx_ovl_dg + 49);

    auto ts_yy_xxzz = pbuffer.data(idx_ovl_dg + 50);

    auto ts_yy_xyyy = pbuffer.data(idx_ovl_dg + 51);

    auto ts_yy_xyyz = pbuffer.data(idx_ovl_dg + 52);

    auto ts_yy_xyzz = pbuffer.data(idx_ovl_dg + 53);

    auto ts_yy_xzzz = pbuffer.data(idx_ovl_dg + 54);

    auto ts_yy_yyyy = pbuffer.data(idx_ovl_dg + 55);

    auto ts_yy_yyyz = pbuffer.data(idx_ovl_dg + 56);

    auto ts_yy_yyzz = pbuffer.data(idx_ovl_dg + 57);

    auto ts_yy_yzzz = pbuffer.data(idx_ovl_dg + 58);

    auto ts_yy_zzzz = pbuffer.data(idx_ovl_dg + 59);

    auto ts_yz_xxxz = pbuffer.data(idx_ovl_dg + 62);

    auto ts_yz_xxyz = pbuffer.data(idx_ovl_dg + 64);

    auto ts_yz_xxzz = pbuffer.data(idx_ovl_dg + 65);

    auto ts_yz_xyyz = pbuffer.data(idx_ovl_dg + 67);

    auto ts_yz_xyzz = pbuffer.data(idx_ovl_dg + 68);

    auto ts_yz_xzzz = pbuffer.data(idx_ovl_dg + 69);

    auto ts_yz_yyyy = pbuffer.data(idx_ovl_dg + 70);

    auto ts_yz_yyyz = pbuffer.data(idx_ovl_dg + 71);

    auto ts_yz_yyzz = pbuffer.data(idx_ovl_dg + 72);

    auto ts_yz_yzzz = pbuffer.data(idx_ovl_dg + 73);

    auto ts_yz_zzzz = pbuffer.data(idx_ovl_dg + 74);

    auto ts_zz_xxxx = pbuffer.data(idx_ovl_dg + 75);

    auto ts_zz_xxxy = pbuffer.data(idx_ovl_dg + 76);

    auto ts_zz_xxxz = pbuffer.data(idx_ovl_dg + 77);

    auto ts_zz_xxyy = pbuffer.data(idx_ovl_dg + 78);

    auto ts_zz_xxyz = pbuffer.data(idx_ovl_dg + 79);

    auto ts_zz_xxzz = pbuffer.data(idx_ovl_dg + 80);

    auto ts_zz_xyyy = pbuffer.data(idx_ovl_dg + 81);

    auto ts_zz_xyyz = pbuffer.data(idx_ovl_dg + 82);

    auto ts_zz_xyzz = pbuffer.data(idx_ovl_dg + 83);

    auto ts_zz_xzzz = pbuffer.data(idx_ovl_dg + 84);

    auto ts_zz_yyyy = pbuffer.data(idx_ovl_dg + 85);

    auto ts_zz_yyyz = pbuffer.data(idx_ovl_dg + 86);

    auto ts_zz_yyzz = pbuffer.data(idx_ovl_dg + 87);

    auto ts_zz_yzzz = pbuffer.data(idx_ovl_dg + 88);

    auto ts_zz_zzzz = pbuffer.data(idx_ovl_dg + 89);

    // Set up 0-15 components of targeted buffer : FG

    auto ts_xxx_xxxx = pbuffer.data(idx_ovl_fg);

    auto ts_xxx_xxxy = pbuffer.data(idx_ovl_fg + 1);

    auto ts_xxx_xxxz = pbuffer.data(idx_ovl_fg + 2);

    auto ts_xxx_xxyy = pbuffer.data(idx_ovl_fg + 3);

    auto ts_xxx_xxyz = pbuffer.data(idx_ovl_fg + 4);

    auto ts_xxx_xxzz = pbuffer.data(idx_ovl_fg + 5);

    auto ts_xxx_xyyy = pbuffer.data(idx_ovl_fg + 6);

    auto ts_xxx_xyyz = pbuffer.data(idx_ovl_fg + 7);

    auto ts_xxx_xyzz = pbuffer.data(idx_ovl_fg + 8);

    auto ts_xxx_xzzz = pbuffer.data(idx_ovl_fg + 9);

    auto ts_xxx_yyyy = pbuffer.data(idx_ovl_fg + 10);

    auto ts_xxx_yyyz = pbuffer.data(idx_ovl_fg + 11);

    auto ts_xxx_yyzz = pbuffer.data(idx_ovl_fg + 12);

    auto ts_xxx_yzzz = pbuffer.data(idx_ovl_fg + 13);

    auto ts_xxx_zzzz = pbuffer.data(idx_ovl_fg + 14);

#pragma omp simd aligned(pa_x,            \
                             ts_x_xxxx,   \
                             ts_x_xxxy,   \
                             ts_x_xxxz,   \
                             ts_x_xxyy,   \
                             ts_x_xxyz,   \
                             ts_x_xxzz,   \
                             ts_x_xyyy,   \
                             ts_x_xyyz,   \
                             ts_x_xyzz,   \
                             ts_x_xzzz,   \
                             ts_x_yyyy,   \
                             ts_x_yyyz,   \
                             ts_x_yyzz,   \
                             ts_x_yzzz,   \
                             ts_x_zzzz,   \
                             ts_xx_xxx,   \
                             ts_xx_xxxx,  \
                             ts_xx_xxxy,  \
                             ts_xx_xxxz,  \
                             ts_xx_xxy,   \
                             ts_xx_xxyy,  \
                             ts_xx_xxyz,  \
                             ts_xx_xxz,   \
                             ts_xx_xxzz,  \
                             ts_xx_xyy,   \
                             ts_xx_xyyy,  \
                             ts_xx_xyyz,  \
                             ts_xx_xyz,   \
                             ts_xx_xyzz,  \
                             ts_xx_xzz,   \
                             ts_xx_xzzz,  \
                             ts_xx_yyy,   \
                             ts_xx_yyyy,  \
                             ts_xx_yyyz,  \
                             ts_xx_yyz,   \
                             ts_xx_yyzz,  \
                             ts_xx_yzz,   \
                             ts_xx_yzzz,  \
                             ts_xx_zzz,   \
                             ts_xx_zzzz,  \
                             ts_xxx_xxxx, \
                             ts_xxx_xxxy, \
                             ts_xxx_xxxz, \
                             ts_xxx_xxyy, \
                             ts_xxx_xxyz, \
                             ts_xxx_xxzz, \
                             ts_xxx_xyyy, \
                             ts_xxx_xyyz, \
                             ts_xxx_xyzz, \
                             ts_xxx_xzzz, \
                             ts_xxx_yyyy, \
                             ts_xxx_yyyz, \
                             ts_xxx_yyzz, \
                             ts_xxx_yzzz, \
                             ts_xxx_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxx_xxxx[i] = 2.0 * ts_x_xxxx[i] * fe_0 + 4.0 * ts_xx_xxx[i] * fe_0 + ts_xx_xxxx[i] * pa_x[i];

        ts_xxx_xxxy[i] = 2.0 * ts_x_xxxy[i] * fe_0 + 3.0 * ts_xx_xxy[i] * fe_0 + ts_xx_xxxy[i] * pa_x[i];

        ts_xxx_xxxz[i] = 2.0 * ts_x_xxxz[i] * fe_0 + 3.0 * ts_xx_xxz[i] * fe_0 + ts_xx_xxxz[i] * pa_x[i];

        ts_xxx_xxyy[i] = 2.0 * ts_x_xxyy[i] * fe_0 + 2.0 * ts_xx_xyy[i] * fe_0 + ts_xx_xxyy[i] * pa_x[i];

        ts_xxx_xxyz[i] = 2.0 * ts_x_xxyz[i] * fe_0 + 2.0 * ts_xx_xyz[i] * fe_0 + ts_xx_xxyz[i] * pa_x[i];

        ts_xxx_xxzz[i] = 2.0 * ts_x_xxzz[i] * fe_0 + 2.0 * ts_xx_xzz[i] * fe_0 + ts_xx_xxzz[i] * pa_x[i];

        ts_xxx_xyyy[i] = 2.0 * ts_x_xyyy[i] * fe_0 + ts_xx_yyy[i] * fe_0 + ts_xx_xyyy[i] * pa_x[i];

        ts_xxx_xyyz[i] = 2.0 * ts_x_xyyz[i] * fe_0 + ts_xx_yyz[i] * fe_0 + ts_xx_xyyz[i] * pa_x[i];

        ts_xxx_xyzz[i] = 2.0 * ts_x_xyzz[i] * fe_0 + ts_xx_yzz[i] * fe_0 + ts_xx_xyzz[i] * pa_x[i];

        ts_xxx_xzzz[i] = 2.0 * ts_x_xzzz[i] * fe_0 + ts_xx_zzz[i] * fe_0 + ts_xx_xzzz[i] * pa_x[i];

        ts_xxx_yyyy[i] = 2.0 * ts_x_yyyy[i] * fe_0 + ts_xx_yyyy[i] * pa_x[i];

        ts_xxx_yyyz[i] = 2.0 * ts_x_yyyz[i] * fe_0 + ts_xx_yyyz[i] * pa_x[i];

        ts_xxx_yyzz[i] = 2.0 * ts_x_yyzz[i] * fe_0 + ts_xx_yyzz[i] * pa_x[i];

        ts_xxx_yzzz[i] = 2.0 * ts_x_yzzz[i] * fe_0 + ts_xx_yzzz[i] * pa_x[i];

        ts_xxx_zzzz[i] = 2.0 * ts_x_zzzz[i] * fe_0 + ts_xx_zzzz[i] * pa_x[i];
    }

    // Set up 15-30 components of targeted buffer : FG

    auto ts_xxy_xxxx = pbuffer.data(idx_ovl_fg + 15);

    auto ts_xxy_xxxy = pbuffer.data(idx_ovl_fg + 16);

    auto ts_xxy_xxxz = pbuffer.data(idx_ovl_fg + 17);

    auto ts_xxy_xxyy = pbuffer.data(idx_ovl_fg + 18);

    auto ts_xxy_xxyz = pbuffer.data(idx_ovl_fg + 19);

    auto ts_xxy_xxzz = pbuffer.data(idx_ovl_fg + 20);

    auto ts_xxy_xyyy = pbuffer.data(idx_ovl_fg + 21);

    auto ts_xxy_xyyz = pbuffer.data(idx_ovl_fg + 22);

    auto ts_xxy_xyzz = pbuffer.data(idx_ovl_fg + 23);

    auto ts_xxy_xzzz = pbuffer.data(idx_ovl_fg + 24);

    auto ts_xxy_yyyy = pbuffer.data(idx_ovl_fg + 25);

    auto ts_xxy_yyyz = pbuffer.data(idx_ovl_fg + 26);

    auto ts_xxy_yyzz = pbuffer.data(idx_ovl_fg + 27);

    auto ts_xxy_yzzz = pbuffer.data(idx_ovl_fg + 28);

    auto ts_xxy_zzzz = pbuffer.data(idx_ovl_fg + 29);

#pragma omp simd aligned(pa_x,            \
                             pa_y,        \
                             ts_xx_xxx,   \
                             ts_xx_xxxx,  \
                             ts_xx_xxxy,  \
                             ts_xx_xxxz,  \
                             ts_xx_xxy,   \
                             ts_xx_xxyy,  \
                             ts_xx_xxyz,  \
                             ts_xx_xxz,   \
                             ts_xx_xxzz,  \
                             ts_xx_xyy,   \
                             ts_xx_xyyy,  \
                             ts_xx_xyyz,  \
                             ts_xx_xyz,   \
                             ts_xx_xyzz,  \
                             ts_xx_xzz,   \
                             ts_xx_xzzz,  \
                             ts_xx_zzzz,  \
                             ts_xxy_xxxx, \
                             ts_xxy_xxxy, \
                             ts_xxy_xxxz, \
                             ts_xxy_xxyy, \
                             ts_xxy_xxyz, \
                             ts_xxy_xxzz, \
                             ts_xxy_xyyy, \
                             ts_xxy_xyyz, \
                             ts_xxy_xyzz, \
                             ts_xxy_xzzz, \
                             ts_xxy_yyyy, \
                             ts_xxy_yyyz, \
                             ts_xxy_yyzz, \
                             ts_xxy_yzzz, \
                             ts_xxy_zzzz, \
                             ts_xy_yyyy,  \
                             ts_xy_yyyz,  \
                             ts_xy_yyzz,  \
                             ts_xy_yzzz,  \
                             ts_y_yyyy,   \
                             ts_y_yyyz,   \
                             ts_y_yyzz,   \
                             ts_y_yzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxy_xxxx[i] = ts_xx_xxxx[i] * pa_y[i];

        ts_xxy_xxxy[i] = ts_xx_xxx[i] * fe_0 + ts_xx_xxxy[i] * pa_y[i];

        ts_xxy_xxxz[i] = ts_xx_xxxz[i] * pa_y[i];

        ts_xxy_xxyy[i] = 2.0 * ts_xx_xxy[i] * fe_0 + ts_xx_xxyy[i] * pa_y[i];

        ts_xxy_xxyz[i] = ts_xx_xxz[i] * fe_0 + ts_xx_xxyz[i] * pa_y[i];

        ts_xxy_xxzz[i] = ts_xx_xxzz[i] * pa_y[i];

        ts_xxy_xyyy[i] = 3.0 * ts_xx_xyy[i] * fe_0 + ts_xx_xyyy[i] * pa_y[i];

        ts_xxy_xyyz[i] = 2.0 * ts_xx_xyz[i] * fe_0 + ts_xx_xyyz[i] * pa_y[i];

        ts_xxy_xyzz[i] = ts_xx_xzz[i] * fe_0 + ts_xx_xyzz[i] * pa_y[i];

        ts_xxy_xzzz[i] = ts_xx_xzzz[i] * pa_y[i];

        ts_xxy_yyyy[i] = ts_y_yyyy[i] * fe_0 + ts_xy_yyyy[i] * pa_x[i];

        ts_xxy_yyyz[i] = ts_y_yyyz[i] * fe_0 + ts_xy_yyyz[i] * pa_x[i];

        ts_xxy_yyzz[i] = ts_y_yyzz[i] * fe_0 + ts_xy_yyzz[i] * pa_x[i];

        ts_xxy_yzzz[i] = ts_y_yzzz[i] * fe_0 + ts_xy_yzzz[i] * pa_x[i];

        ts_xxy_zzzz[i] = ts_xx_zzzz[i] * pa_y[i];
    }

    // Set up 30-45 components of targeted buffer : FG

    auto ts_xxz_xxxx = pbuffer.data(idx_ovl_fg + 30);

    auto ts_xxz_xxxy = pbuffer.data(idx_ovl_fg + 31);

    auto ts_xxz_xxxz = pbuffer.data(idx_ovl_fg + 32);

    auto ts_xxz_xxyy = pbuffer.data(idx_ovl_fg + 33);

    auto ts_xxz_xxyz = pbuffer.data(idx_ovl_fg + 34);

    auto ts_xxz_xxzz = pbuffer.data(idx_ovl_fg + 35);

    auto ts_xxz_xyyy = pbuffer.data(idx_ovl_fg + 36);

    auto ts_xxz_xyyz = pbuffer.data(idx_ovl_fg + 37);

    auto ts_xxz_xyzz = pbuffer.data(idx_ovl_fg + 38);

    auto ts_xxz_xzzz = pbuffer.data(idx_ovl_fg + 39);

    auto ts_xxz_yyyy = pbuffer.data(idx_ovl_fg + 40);

    auto ts_xxz_yyyz = pbuffer.data(idx_ovl_fg + 41);

    auto ts_xxz_yyzz = pbuffer.data(idx_ovl_fg + 42);

    auto ts_xxz_yzzz = pbuffer.data(idx_ovl_fg + 43);

    auto ts_xxz_zzzz = pbuffer.data(idx_ovl_fg + 44);

#pragma omp simd aligned(pa_x,            \
                             pa_z,        \
                             ts_xx_xxx,   \
                             ts_xx_xxxx,  \
                             ts_xx_xxxy,  \
                             ts_xx_xxxz,  \
                             ts_xx_xxy,   \
                             ts_xx_xxyy,  \
                             ts_xx_xxyz,  \
                             ts_xx_xxz,   \
                             ts_xx_xxzz,  \
                             ts_xx_xyy,   \
                             ts_xx_xyyy,  \
                             ts_xx_xyyz,  \
                             ts_xx_xyz,   \
                             ts_xx_xyzz,  \
                             ts_xx_xzz,   \
                             ts_xx_xzzz,  \
                             ts_xx_yyyy,  \
                             ts_xxz_xxxx, \
                             ts_xxz_xxxy, \
                             ts_xxz_xxxz, \
                             ts_xxz_xxyy, \
                             ts_xxz_xxyz, \
                             ts_xxz_xxzz, \
                             ts_xxz_xyyy, \
                             ts_xxz_xyyz, \
                             ts_xxz_xyzz, \
                             ts_xxz_xzzz, \
                             ts_xxz_yyyy, \
                             ts_xxz_yyyz, \
                             ts_xxz_yyzz, \
                             ts_xxz_yzzz, \
                             ts_xxz_zzzz, \
                             ts_xz_yyyz,  \
                             ts_xz_yyzz,  \
                             ts_xz_yzzz,  \
                             ts_xz_zzzz,  \
                             ts_z_yyyz,   \
                             ts_z_yyzz,   \
                             ts_z_yzzz,   \
                             ts_z_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxz_xxxx[i] = ts_xx_xxxx[i] * pa_z[i];

        ts_xxz_xxxy[i] = ts_xx_xxxy[i] * pa_z[i];

        ts_xxz_xxxz[i] = ts_xx_xxx[i] * fe_0 + ts_xx_xxxz[i] * pa_z[i];

        ts_xxz_xxyy[i] = ts_xx_xxyy[i] * pa_z[i];

        ts_xxz_xxyz[i] = ts_xx_xxy[i] * fe_0 + ts_xx_xxyz[i] * pa_z[i];

        ts_xxz_xxzz[i] = 2.0 * ts_xx_xxz[i] * fe_0 + ts_xx_xxzz[i] * pa_z[i];

        ts_xxz_xyyy[i] = ts_xx_xyyy[i] * pa_z[i];

        ts_xxz_xyyz[i] = ts_xx_xyy[i] * fe_0 + ts_xx_xyyz[i] * pa_z[i];

        ts_xxz_xyzz[i] = 2.0 * ts_xx_xyz[i] * fe_0 + ts_xx_xyzz[i] * pa_z[i];

        ts_xxz_xzzz[i] = 3.0 * ts_xx_xzz[i] * fe_0 + ts_xx_xzzz[i] * pa_z[i];

        ts_xxz_yyyy[i] = ts_xx_yyyy[i] * pa_z[i];

        ts_xxz_yyyz[i] = ts_z_yyyz[i] * fe_0 + ts_xz_yyyz[i] * pa_x[i];

        ts_xxz_yyzz[i] = ts_z_yyzz[i] * fe_0 + ts_xz_yyzz[i] * pa_x[i];

        ts_xxz_yzzz[i] = ts_z_yzzz[i] * fe_0 + ts_xz_yzzz[i] * pa_x[i];

        ts_xxz_zzzz[i] = ts_z_zzzz[i] * fe_0 + ts_xz_zzzz[i] * pa_x[i];
    }

    // Set up 45-60 components of targeted buffer : FG

    auto ts_xyy_xxxx = pbuffer.data(idx_ovl_fg + 45);

    auto ts_xyy_xxxy = pbuffer.data(idx_ovl_fg + 46);

    auto ts_xyy_xxxz = pbuffer.data(idx_ovl_fg + 47);

    auto ts_xyy_xxyy = pbuffer.data(idx_ovl_fg + 48);

    auto ts_xyy_xxyz = pbuffer.data(idx_ovl_fg + 49);

    auto ts_xyy_xxzz = pbuffer.data(idx_ovl_fg + 50);

    auto ts_xyy_xyyy = pbuffer.data(idx_ovl_fg + 51);

    auto ts_xyy_xyyz = pbuffer.data(idx_ovl_fg + 52);

    auto ts_xyy_xyzz = pbuffer.data(idx_ovl_fg + 53);

    auto ts_xyy_xzzz = pbuffer.data(idx_ovl_fg + 54);

    auto ts_xyy_yyyy = pbuffer.data(idx_ovl_fg + 55);

    auto ts_xyy_yyyz = pbuffer.data(idx_ovl_fg + 56);

    auto ts_xyy_yyzz = pbuffer.data(idx_ovl_fg + 57);

    auto ts_xyy_yzzz = pbuffer.data(idx_ovl_fg + 58);

    auto ts_xyy_zzzz = pbuffer.data(idx_ovl_fg + 59);

#pragma omp simd aligned(pa_x,            \
                             ts_xyy_xxxx, \
                             ts_xyy_xxxy, \
                             ts_xyy_xxxz, \
                             ts_xyy_xxyy, \
                             ts_xyy_xxyz, \
                             ts_xyy_xxzz, \
                             ts_xyy_xyyy, \
                             ts_xyy_xyyz, \
                             ts_xyy_xyzz, \
                             ts_xyy_xzzz, \
                             ts_xyy_yyyy, \
                             ts_xyy_yyyz, \
                             ts_xyy_yyzz, \
                             ts_xyy_yzzz, \
                             ts_xyy_zzzz, \
                             ts_yy_xxx,   \
                             ts_yy_xxxx,  \
                             ts_yy_xxxy,  \
                             ts_yy_xxxz,  \
                             ts_yy_xxy,   \
                             ts_yy_xxyy,  \
                             ts_yy_xxyz,  \
                             ts_yy_xxz,   \
                             ts_yy_xxzz,  \
                             ts_yy_xyy,   \
                             ts_yy_xyyy,  \
                             ts_yy_xyyz,  \
                             ts_yy_xyz,   \
                             ts_yy_xyzz,  \
                             ts_yy_xzz,   \
                             ts_yy_xzzz,  \
                             ts_yy_yyy,   \
                             ts_yy_yyyy,  \
                             ts_yy_yyyz,  \
                             ts_yy_yyz,   \
                             ts_yy_yyzz,  \
                             ts_yy_yzz,   \
                             ts_yy_yzzz,  \
                             ts_yy_zzz,   \
                             ts_yy_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyy_xxxx[i] = 4.0 * ts_yy_xxx[i] * fe_0 + ts_yy_xxxx[i] * pa_x[i];

        ts_xyy_xxxy[i] = 3.0 * ts_yy_xxy[i] * fe_0 + ts_yy_xxxy[i] * pa_x[i];

        ts_xyy_xxxz[i] = 3.0 * ts_yy_xxz[i] * fe_0 + ts_yy_xxxz[i] * pa_x[i];

        ts_xyy_xxyy[i] = 2.0 * ts_yy_xyy[i] * fe_0 + ts_yy_xxyy[i] * pa_x[i];

        ts_xyy_xxyz[i] = 2.0 * ts_yy_xyz[i] * fe_0 + ts_yy_xxyz[i] * pa_x[i];

        ts_xyy_xxzz[i] = 2.0 * ts_yy_xzz[i] * fe_0 + ts_yy_xxzz[i] * pa_x[i];

        ts_xyy_xyyy[i] = ts_yy_yyy[i] * fe_0 + ts_yy_xyyy[i] * pa_x[i];

        ts_xyy_xyyz[i] = ts_yy_yyz[i] * fe_0 + ts_yy_xyyz[i] * pa_x[i];

        ts_xyy_xyzz[i] = ts_yy_yzz[i] * fe_0 + ts_yy_xyzz[i] * pa_x[i];

        ts_xyy_xzzz[i] = ts_yy_zzz[i] * fe_0 + ts_yy_xzzz[i] * pa_x[i];

        ts_xyy_yyyy[i] = ts_yy_yyyy[i] * pa_x[i];

        ts_xyy_yyyz[i] = ts_yy_yyyz[i] * pa_x[i];

        ts_xyy_yyzz[i] = ts_yy_yyzz[i] * pa_x[i];

        ts_xyy_yzzz[i] = ts_yy_yzzz[i] * pa_x[i];

        ts_xyy_zzzz[i] = ts_yy_zzzz[i] * pa_x[i];
    }

    // Set up 60-75 components of targeted buffer : FG

    auto ts_xyz_xxxx = pbuffer.data(idx_ovl_fg + 60);

    auto ts_xyz_xxxy = pbuffer.data(idx_ovl_fg + 61);

    auto ts_xyz_xxxz = pbuffer.data(idx_ovl_fg + 62);

    auto ts_xyz_xxyy = pbuffer.data(idx_ovl_fg + 63);

    auto ts_xyz_xxyz = pbuffer.data(idx_ovl_fg + 64);

    auto ts_xyz_xxzz = pbuffer.data(idx_ovl_fg + 65);

    auto ts_xyz_xyyy = pbuffer.data(idx_ovl_fg + 66);

    auto ts_xyz_xyyz = pbuffer.data(idx_ovl_fg + 67);

    auto ts_xyz_xyzz = pbuffer.data(idx_ovl_fg + 68);

    auto ts_xyz_xzzz = pbuffer.data(idx_ovl_fg + 69);

    auto ts_xyz_yyyy = pbuffer.data(idx_ovl_fg + 70);

    auto ts_xyz_yyyz = pbuffer.data(idx_ovl_fg + 71);

    auto ts_xyz_yyzz = pbuffer.data(idx_ovl_fg + 72);

    auto ts_xyz_yzzz = pbuffer.data(idx_ovl_fg + 73);

    auto ts_xyz_zzzz = pbuffer.data(idx_ovl_fg + 74);

#pragma omp simd aligned(pa_x,            \
                             pa_y,        \
                             pa_z,        \
                             ts_xy_xxxy,  \
                             ts_xy_xxyy,  \
                             ts_xy_xyyy,  \
                             ts_xyz_xxxx, \
                             ts_xyz_xxxy, \
                             ts_xyz_xxxz, \
                             ts_xyz_xxyy, \
                             ts_xyz_xxyz, \
                             ts_xyz_xxzz, \
                             ts_xyz_xyyy, \
                             ts_xyz_xyyz, \
                             ts_xyz_xyzz, \
                             ts_xyz_xzzz, \
                             ts_xyz_yyyy, \
                             ts_xyz_yyyz, \
                             ts_xyz_yyzz, \
                             ts_xyz_yzzz, \
                             ts_xyz_zzzz, \
                             ts_xz_xxxx,  \
                             ts_xz_xxxz,  \
                             ts_xz_xxzz,  \
                             ts_xz_xzzz,  \
                             ts_yz_xxyz,  \
                             ts_yz_xyyz,  \
                             ts_yz_xyz,   \
                             ts_yz_xyzz,  \
                             ts_yz_yyyy,  \
                             ts_yz_yyyz,  \
                             ts_yz_yyz,   \
                             ts_yz_yyzz,  \
                             ts_yz_yzz,   \
                             ts_yz_yzzz,  \
                             ts_yz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyz_xxxx[i] = ts_xz_xxxx[i] * pa_y[i];

        ts_xyz_xxxy[i] = ts_xy_xxxy[i] * pa_z[i];

        ts_xyz_xxxz[i] = ts_xz_xxxz[i] * pa_y[i];

        ts_xyz_xxyy[i] = ts_xy_xxyy[i] * pa_z[i];

        ts_xyz_xxyz[i] = 2.0 * ts_yz_xyz[i] * fe_0 + ts_yz_xxyz[i] * pa_x[i];

        ts_xyz_xxzz[i] = ts_xz_xxzz[i] * pa_y[i];

        ts_xyz_xyyy[i] = ts_xy_xyyy[i] * pa_z[i];

        ts_xyz_xyyz[i] = ts_yz_yyz[i] * fe_0 + ts_yz_xyyz[i] * pa_x[i];

        ts_xyz_xyzz[i] = ts_yz_yzz[i] * fe_0 + ts_yz_xyzz[i] * pa_x[i];

        ts_xyz_xzzz[i] = ts_xz_xzzz[i] * pa_y[i];

        ts_xyz_yyyy[i] = ts_yz_yyyy[i] * pa_x[i];

        ts_xyz_yyyz[i] = ts_yz_yyyz[i] * pa_x[i];

        ts_xyz_yyzz[i] = ts_yz_yyzz[i] * pa_x[i];

        ts_xyz_yzzz[i] = ts_yz_yzzz[i] * pa_x[i];

        ts_xyz_zzzz[i] = ts_yz_zzzz[i] * pa_x[i];
    }

    // Set up 75-90 components of targeted buffer : FG

    auto ts_xzz_xxxx = pbuffer.data(idx_ovl_fg + 75);

    auto ts_xzz_xxxy = pbuffer.data(idx_ovl_fg + 76);

    auto ts_xzz_xxxz = pbuffer.data(idx_ovl_fg + 77);

    auto ts_xzz_xxyy = pbuffer.data(idx_ovl_fg + 78);

    auto ts_xzz_xxyz = pbuffer.data(idx_ovl_fg + 79);

    auto ts_xzz_xxzz = pbuffer.data(idx_ovl_fg + 80);

    auto ts_xzz_xyyy = pbuffer.data(idx_ovl_fg + 81);

    auto ts_xzz_xyyz = pbuffer.data(idx_ovl_fg + 82);

    auto ts_xzz_xyzz = pbuffer.data(idx_ovl_fg + 83);

    auto ts_xzz_xzzz = pbuffer.data(idx_ovl_fg + 84);

    auto ts_xzz_yyyy = pbuffer.data(idx_ovl_fg + 85);

    auto ts_xzz_yyyz = pbuffer.data(idx_ovl_fg + 86);

    auto ts_xzz_yyzz = pbuffer.data(idx_ovl_fg + 87);

    auto ts_xzz_yzzz = pbuffer.data(idx_ovl_fg + 88);

    auto ts_xzz_zzzz = pbuffer.data(idx_ovl_fg + 89);

#pragma omp simd aligned(pa_x,            \
                             ts_xzz_xxxx, \
                             ts_xzz_xxxy, \
                             ts_xzz_xxxz, \
                             ts_xzz_xxyy, \
                             ts_xzz_xxyz, \
                             ts_xzz_xxzz, \
                             ts_xzz_xyyy, \
                             ts_xzz_xyyz, \
                             ts_xzz_xyzz, \
                             ts_xzz_xzzz, \
                             ts_xzz_yyyy, \
                             ts_xzz_yyyz, \
                             ts_xzz_yyzz, \
                             ts_xzz_yzzz, \
                             ts_xzz_zzzz, \
                             ts_zz_xxx,   \
                             ts_zz_xxxx,  \
                             ts_zz_xxxy,  \
                             ts_zz_xxxz,  \
                             ts_zz_xxy,   \
                             ts_zz_xxyy,  \
                             ts_zz_xxyz,  \
                             ts_zz_xxz,   \
                             ts_zz_xxzz,  \
                             ts_zz_xyy,   \
                             ts_zz_xyyy,  \
                             ts_zz_xyyz,  \
                             ts_zz_xyz,   \
                             ts_zz_xyzz,  \
                             ts_zz_xzz,   \
                             ts_zz_xzzz,  \
                             ts_zz_yyy,   \
                             ts_zz_yyyy,  \
                             ts_zz_yyyz,  \
                             ts_zz_yyz,   \
                             ts_zz_yyzz,  \
                             ts_zz_yzz,   \
                             ts_zz_yzzz,  \
                             ts_zz_zzz,   \
                             ts_zz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xzz_xxxx[i] = 4.0 * ts_zz_xxx[i] * fe_0 + ts_zz_xxxx[i] * pa_x[i];

        ts_xzz_xxxy[i] = 3.0 * ts_zz_xxy[i] * fe_0 + ts_zz_xxxy[i] * pa_x[i];

        ts_xzz_xxxz[i] = 3.0 * ts_zz_xxz[i] * fe_0 + ts_zz_xxxz[i] * pa_x[i];

        ts_xzz_xxyy[i] = 2.0 * ts_zz_xyy[i] * fe_0 + ts_zz_xxyy[i] * pa_x[i];

        ts_xzz_xxyz[i] = 2.0 * ts_zz_xyz[i] * fe_0 + ts_zz_xxyz[i] * pa_x[i];

        ts_xzz_xxzz[i] = 2.0 * ts_zz_xzz[i] * fe_0 + ts_zz_xxzz[i] * pa_x[i];

        ts_xzz_xyyy[i] = ts_zz_yyy[i] * fe_0 + ts_zz_xyyy[i] * pa_x[i];

        ts_xzz_xyyz[i] = ts_zz_yyz[i] * fe_0 + ts_zz_xyyz[i] * pa_x[i];

        ts_xzz_xyzz[i] = ts_zz_yzz[i] * fe_0 + ts_zz_xyzz[i] * pa_x[i];

        ts_xzz_xzzz[i] = ts_zz_zzz[i] * fe_0 + ts_zz_xzzz[i] * pa_x[i];

        ts_xzz_yyyy[i] = ts_zz_yyyy[i] * pa_x[i];

        ts_xzz_yyyz[i] = ts_zz_yyyz[i] * pa_x[i];

        ts_xzz_yyzz[i] = ts_zz_yyzz[i] * pa_x[i];

        ts_xzz_yzzz[i] = ts_zz_yzzz[i] * pa_x[i];

        ts_xzz_zzzz[i] = ts_zz_zzzz[i] * pa_x[i];
    }

    // Set up 90-105 components of targeted buffer : FG

    auto ts_yyy_xxxx = pbuffer.data(idx_ovl_fg + 90);

    auto ts_yyy_xxxy = pbuffer.data(idx_ovl_fg + 91);

    auto ts_yyy_xxxz = pbuffer.data(idx_ovl_fg + 92);

    auto ts_yyy_xxyy = pbuffer.data(idx_ovl_fg + 93);

    auto ts_yyy_xxyz = pbuffer.data(idx_ovl_fg + 94);

    auto ts_yyy_xxzz = pbuffer.data(idx_ovl_fg + 95);

    auto ts_yyy_xyyy = pbuffer.data(idx_ovl_fg + 96);

    auto ts_yyy_xyyz = pbuffer.data(idx_ovl_fg + 97);

    auto ts_yyy_xyzz = pbuffer.data(idx_ovl_fg + 98);

    auto ts_yyy_xzzz = pbuffer.data(idx_ovl_fg + 99);

    auto ts_yyy_yyyy = pbuffer.data(idx_ovl_fg + 100);

    auto ts_yyy_yyyz = pbuffer.data(idx_ovl_fg + 101);

    auto ts_yyy_yyzz = pbuffer.data(idx_ovl_fg + 102);

    auto ts_yyy_yzzz = pbuffer.data(idx_ovl_fg + 103);

    auto ts_yyy_zzzz = pbuffer.data(idx_ovl_fg + 104);

#pragma omp simd aligned(pa_y,            \
                             ts_y_xxxx,   \
                             ts_y_xxxy,   \
                             ts_y_xxxz,   \
                             ts_y_xxyy,   \
                             ts_y_xxyz,   \
                             ts_y_xxzz,   \
                             ts_y_xyyy,   \
                             ts_y_xyyz,   \
                             ts_y_xyzz,   \
                             ts_y_xzzz,   \
                             ts_y_yyyy,   \
                             ts_y_yyyz,   \
                             ts_y_yyzz,   \
                             ts_y_yzzz,   \
                             ts_y_zzzz,   \
                             ts_yy_xxx,   \
                             ts_yy_xxxx,  \
                             ts_yy_xxxy,  \
                             ts_yy_xxxz,  \
                             ts_yy_xxy,   \
                             ts_yy_xxyy,  \
                             ts_yy_xxyz,  \
                             ts_yy_xxz,   \
                             ts_yy_xxzz,  \
                             ts_yy_xyy,   \
                             ts_yy_xyyy,  \
                             ts_yy_xyyz,  \
                             ts_yy_xyz,   \
                             ts_yy_xyzz,  \
                             ts_yy_xzz,   \
                             ts_yy_xzzz,  \
                             ts_yy_yyy,   \
                             ts_yy_yyyy,  \
                             ts_yy_yyyz,  \
                             ts_yy_yyz,   \
                             ts_yy_yyzz,  \
                             ts_yy_yzz,   \
                             ts_yy_yzzz,  \
                             ts_yy_zzz,   \
                             ts_yy_zzzz,  \
                             ts_yyy_xxxx, \
                             ts_yyy_xxxy, \
                             ts_yyy_xxxz, \
                             ts_yyy_xxyy, \
                             ts_yyy_xxyz, \
                             ts_yyy_xxzz, \
                             ts_yyy_xyyy, \
                             ts_yyy_xyyz, \
                             ts_yyy_xyzz, \
                             ts_yyy_xzzz, \
                             ts_yyy_yyyy, \
                             ts_yyy_yyyz, \
                             ts_yyy_yyzz, \
                             ts_yyy_yzzz, \
                             ts_yyy_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyy_xxxx[i] = 2.0 * ts_y_xxxx[i] * fe_0 + ts_yy_xxxx[i] * pa_y[i];

        ts_yyy_xxxy[i] = 2.0 * ts_y_xxxy[i] * fe_0 + ts_yy_xxx[i] * fe_0 + ts_yy_xxxy[i] * pa_y[i];

        ts_yyy_xxxz[i] = 2.0 * ts_y_xxxz[i] * fe_0 + ts_yy_xxxz[i] * pa_y[i];

        ts_yyy_xxyy[i] = 2.0 * ts_y_xxyy[i] * fe_0 + 2.0 * ts_yy_xxy[i] * fe_0 + ts_yy_xxyy[i] * pa_y[i];

        ts_yyy_xxyz[i] = 2.0 * ts_y_xxyz[i] * fe_0 + ts_yy_xxz[i] * fe_0 + ts_yy_xxyz[i] * pa_y[i];

        ts_yyy_xxzz[i] = 2.0 * ts_y_xxzz[i] * fe_0 + ts_yy_xxzz[i] * pa_y[i];

        ts_yyy_xyyy[i] = 2.0 * ts_y_xyyy[i] * fe_0 + 3.0 * ts_yy_xyy[i] * fe_0 + ts_yy_xyyy[i] * pa_y[i];

        ts_yyy_xyyz[i] = 2.0 * ts_y_xyyz[i] * fe_0 + 2.0 * ts_yy_xyz[i] * fe_0 + ts_yy_xyyz[i] * pa_y[i];

        ts_yyy_xyzz[i] = 2.0 * ts_y_xyzz[i] * fe_0 + ts_yy_xzz[i] * fe_0 + ts_yy_xyzz[i] * pa_y[i];

        ts_yyy_xzzz[i] = 2.0 * ts_y_xzzz[i] * fe_0 + ts_yy_xzzz[i] * pa_y[i];

        ts_yyy_yyyy[i] = 2.0 * ts_y_yyyy[i] * fe_0 + 4.0 * ts_yy_yyy[i] * fe_0 + ts_yy_yyyy[i] * pa_y[i];

        ts_yyy_yyyz[i] = 2.0 * ts_y_yyyz[i] * fe_0 + 3.0 * ts_yy_yyz[i] * fe_0 + ts_yy_yyyz[i] * pa_y[i];

        ts_yyy_yyzz[i] = 2.0 * ts_y_yyzz[i] * fe_0 + 2.0 * ts_yy_yzz[i] * fe_0 + ts_yy_yyzz[i] * pa_y[i];

        ts_yyy_yzzz[i] = 2.0 * ts_y_yzzz[i] * fe_0 + ts_yy_zzz[i] * fe_0 + ts_yy_yzzz[i] * pa_y[i];

        ts_yyy_zzzz[i] = 2.0 * ts_y_zzzz[i] * fe_0 + ts_yy_zzzz[i] * pa_y[i];
    }

    // Set up 105-120 components of targeted buffer : FG

    auto ts_yyz_xxxx = pbuffer.data(idx_ovl_fg + 105);

    auto ts_yyz_xxxy = pbuffer.data(idx_ovl_fg + 106);

    auto ts_yyz_xxxz = pbuffer.data(idx_ovl_fg + 107);

    auto ts_yyz_xxyy = pbuffer.data(idx_ovl_fg + 108);

    auto ts_yyz_xxyz = pbuffer.data(idx_ovl_fg + 109);

    auto ts_yyz_xxzz = pbuffer.data(idx_ovl_fg + 110);

    auto ts_yyz_xyyy = pbuffer.data(idx_ovl_fg + 111);

    auto ts_yyz_xyyz = pbuffer.data(idx_ovl_fg + 112);

    auto ts_yyz_xyzz = pbuffer.data(idx_ovl_fg + 113);

    auto ts_yyz_xzzz = pbuffer.data(idx_ovl_fg + 114);

    auto ts_yyz_yyyy = pbuffer.data(idx_ovl_fg + 115);

    auto ts_yyz_yyyz = pbuffer.data(idx_ovl_fg + 116);

    auto ts_yyz_yyzz = pbuffer.data(idx_ovl_fg + 117);

    auto ts_yyz_yzzz = pbuffer.data(idx_ovl_fg + 118);

    auto ts_yyz_zzzz = pbuffer.data(idx_ovl_fg + 119);

#pragma omp simd aligned(pa_y,            \
                             pa_z,        \
                             ts_yy_xxxx,  \
                             ts_yy_xxxy,  \
                             ts_yy_xxy,   \
                             ts_yy_xxyy,  \
                             ts_yy_xxyz,  \
                             ts_yy_xyy,   \
                             ts_yy_xyyy,  \
                             ts_yy_xyyz,  \
                             ts_yy_xyz,   \
                             ts_yy_xyzz,  \
                             ts_yy_yyy,   \
                             ts_yy_yyyy,  \
                             ts_yy_yyyz,  \
                             ts_yy_yyz,   \
                             ts_yy_yyzz,  \
                             ts_yy_yzz,   \
                             ts_yy_yzzz,  \
                             ts_yyz_xxxx, \
                             ts_yyz_xxxy, \
                             ts_yyz_xxxz, \
                             ts_yyz_xxyy, \
                             ts_yyz_xxyz, \
                             ts_yyz_xxzz, \
                             ts_yyz_xyyy, \
                             ts_yyz_xyyz, \
                             ts_yyz_xyzz, \
                             ts_yyz_xzzz, \
                             ts_yyz_yyyy, \
                             ts_yyz_yyyz, \
                             ts_yyz_yyzz, \
                             ts_yyz_yzzz, \
                             ts_yyz_zzzz, \
                             ts_yz_xxxz,  \
                             ts_yz_xxzz,  \
                             ts_yz_xzzz,  \
                             ts_yz_zzzz,  \
                             ts_z_xxxz,   \
                             ts_z_xxzz,   \
                             ts_z_xzzz,   \
                             ts_z_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyz_xxxx[i] = ts_yy_xxxx[i] * pa_z[i];

        ts_yyz_xxxy[i] = ts_yy_xxxy[i] * pa_z[i];

        ts_yyz_xxxz[i] = ts_z_xxxz[i] * fe_0 + ts_yz_xxxz[i] * pa_y[i];

        ts_yyz_xxyy[i] = ts_yy_xxyy[i] * pa_z[i];

        ts_yyz_xxyz[i] = ts_yy_xxy[i] * fe_0 + ts_yy_xxyz[i] * pa_z[i];

        ts_yyz_xxzz[i] = ts_z_xxzz[i] * fe_0 + ts_yz_xxzz[i] * pa_y[i];

        ts_yyz_xyyy[i] = ts_yy_xyyy[i] * pa_z[i];

        ts_yyz_xyyz[i] = ts_yy_xyy[i] * fe_0 + ts_yy_xyyz[i] * pa_z[i];

        ts_yyz_xyzz[i] = 2.0 * ts_yy_xyz[i] * fe_0 + ts_yy_xyzz[i] * pa_z[i];

        ts_yyz_xzzz[i] = ts_z_xzzz[i] * fe_0 + ts_yz_xzzz[i] * pa_y[i];

        ts_yyz_yyyy[i] = ts_yy_yyyy[i] * pa_z[i];

        ts_yyz_yyyz[i] = ts_yy_yyy[i] * fe_0 + ts_yy_yyyz[i] * pa_z[i];

        ts_yyz_yyzz[i] = 2.0 * ts_yy_yyz[i] * fe_0 + ts_yy_yyzz[i] * pa_z[i];

        ts_yyz_yzzz[i] = 3.0 * ts_yy_yzz[i] * fe_0 + ts_yy_yzzz[i] * pa_z[i];

        ts_yyz_zzzz[i] = ts_z_zzzz[i] * fe_0 + ts_yz_zzzz[i] * pa_y[i];
    }

    // Set up 120-135 components of targeted buffer : FG

    auto ts_yzz_xxxx = pbuffer.data(idx_ovl_fg + 120);

    auto ts_yzz_xxxy = pbuffer.data(idx_ovl_fg + 121);

    auto ts_yzz_xxxz = pbuffer.data(idx_ovl_fg + 122);

    auto ts_yzz_xxyy = pbuffer.data(idx_ovl_fg + 123);

    auto ts_yzz_xxyz = pbuffer.data(idx_ovl_fg + 124);

    auto ts_yzz_xxzz = pbuffer.data(idx_ovl_fg + 125);

    auto ts_yzz_xyyy = pbuffer.data(idx_ovl_fg + 126);

    auto ts_yzz_xyyz = pbuffer.data(idx_ovl_fg + 127);

    auto ts_yzz_xyzz = pbuffer.data(idx_ovl_fg + 128);

    auto ts_yzz_xzzz = pbuffer.data(idx_ovl_fg + 129);

    auto ts_yzz_yyyy = pbuffer.data(idx_ovl_fg + 130);

    auto ts_yzz_yyyz = pbuffer.data(idx_ovl_fg + 131);

    auto ts_yzz_yyzz = pbuffer.data(idx_ovl_fg + 132);

    auto ts_yzz_yzzz = pbuffer.data(idx_ovl_fg + 133);

    auto ts_yzz_zzzz = pbuffer.data(idx_ovl_fg + 134);

#pragma omp simd aligned(pa_y,            \
                             ts_yzz_xxxx, \
                             ts_yzz_xxxy, \
                             ts_yzz_xxxz, \
                             ts_yzz_xxyy, \
                             ts_yzz_xxyz, \
                             ts_yzz_xxzz, \
                             ts_yzz_xyyy, \
                             ts_yzz_xyyz, \
                             ts_yzz_xyzz, \
                             ts_yzz_xzzz, \
                             ts_yzz_yyyy, \
                             ts_yzz_yyyz, \
                             ts_yzz_yyzz, \
                             ts_yzz_yzzz, \
                             ts_yzz_zzzz, \
                             ts_zz_xxx,   \
                             ts_zz_xxxx,  \
                             ts_zz_xxxy,  \
                             ts_zz_xxxz,  \
                             ts_zz_xxy,   \
                             ts_zz_xxyy,  \
                             ts_zz_xxyz,  \
                             ts_zz_xxz,   \
                             ts_zz_xxzz,  \
                             ts_zz_xyy,   \
                             ts_zz_xyyy,  \
                             ts_zz_xyyz,  \
                             ts_zz_xyz,   \
                             ts_zz_xyzz,  \
                             ts_zz_xzz,   \
                             ts_zz_xzzz,  \
                             ts_zz_yyy,   \
                             ts_zz_yyyy,  \
                             ts_zz_yyyz,  \
                             ts_zz_yyz,   \
                             ts_zz_yyzz,  \
                             ts_zz_yzz,   \
                             ts_zz_yzzz,  \
                             ts_zz_zzz,   \
                             ts_zz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yzz_xxxx[i] = ts_zz_xxxx[i] * pa_y[i];

        ts_yzz_xxxy[i] = ts_zz_xxx[i] * fe_0 + ts_zz_xxxy[i] * pa_y[i];

        ts_yzz_xxxz[i] = ts_zz_xxxz[i] * pa_y[i];

        ts_yzz_xxyy[i] = 2.0 * ts_zz_xxy[i] * fe_0 + ts_zz_xxyy[i] * pa_y[i];

        ts_yzz_xxyz[i] = ts_zz_xxz[i] * fe_0 + ts_zz_xxyz[i] * pa_y[i];

        ts_yzz_xxzz[i] = ts_zz_xxzz[i] * pa_y[i];

        ts_yzz_xyyy[i] = 3.0 * ts_zz_xyy[i] * fe_0 + ts_zz_xyyy[i] * pa_y[i];

        ts_yzz_xyyz[i] = 2.0 * ts_zz_xyz[i] * fe_0 + ts_zz_xyyz[i] * pa_y[i];

        ts_yzz_xyzz[i] = ts_zz_xzz[i] * fe_0 + ts_zz_xyzz[i] * pa_y[i];

        ts_yzz_xzzz[i] = ts_zz_xzzz[i] * pa_y[i];

        ts_yzz_yyyy[i] = 4.0 * ts_zz_yyy[i] * fe_0 + ts_zz_yyyy[i] * pa_y[i];

        ts_yzz_yyyz[i] = 3.0 * ts_zz_yyz[i] * fe_0 + ts_zz_yyyz[i] * pa_y[i];

        ts_yzz_yyzz[i] = 2.0 * ts_zz_yzz[i] * fe_0 + ts_zz_yyzz[i] * pa_y[i];

        ts_yzz_yzzz[i] = ts_zz_zzz[i] * fe_0 + ts_zz_yzzz[i] * pa_y[i];

        ts_yzz_zzzz[i] = ts_zz_zzzz[i] * pa_y[i];
    }

    // Set up 135-150 components of targeted buffer : FG

    auto ts_zzz_xxxx = pbuffer.data(idx_ovl_fg + 135);

    auto ts_zzz_xxxy = pbuffer.data(idx_ovl_fg + 136);

    auto ts_zzz_xxxz = pbuffer.data(idx_ovl_fg + 137);

    auto ts_zzz_xxyy = pbuffer.data(idx_ovl_fg + 138);

    auto ts_zzz_xxyz = pbuffer.data(idx_ovl_fg + 139);

    auto ts_zzz_xxzz = pbuffer.data(idx_ovl_fg + 140);

    auto ts_zzz_xyyy = pbuffer.data(idx_ovl_fg + 141);

    auto ts_zzz_xyyz = pbuffer.data(idx_ovl_fg + 142);

    auto ts_zzz_xyzz = pbuffer.data(idx_ovl_fg + 143);

    auto ts_zzz_xzzz = pbuffer.data(idx_ovl_fg + 144);

    auto ts_zzz_yyyy = pbuffer.data(idx_ovl_fg + 145);

    auto ts_zzz_yyyz = pbuffer.data(idx_ovl_fg + 146);

    auto ts_zzz_yyzz = pbuffer.data(idx_ovl_fg + 147);

    auto ts_zzz_yzzz = pbuffer.data(idx_ovl_fg + 148);

    auto ts_zzz_zzzz = pbuffer.data(idx_ovl_fg + 149);

#pragma omp simd aligned(pa_z,            \
                             ts_z_xxxx,   \
                             ts_z_xxxy,   \
                             ts_z_xxxz,   \
                             ts_z_xxyy,   \
                             ts_z_xxyz,   \
                             ts_z_xxzz,   \
                             ts_z_xyyy,   \
                             ts_z_xyyz,   \
                             ts_z_xyzz,   \
                             ts_z_xzzz,   \
                             ts_z_yyyy,   \
                             ts_z_yyyz,   \
                             ts_z_yyzz,   \
                             ts_z_yzzz,   \
                             ts_z_zzzz,   \
                             ts_zz_xxx,   \
                             ts_zz_xxxx,  \
                             ts_zz_xxxy,  \
                             ts_zz_xxxz,  \
                             ts_zz_xxy,   \
                             ts_zz_xxyy,  \
                             ts_zz_xxyz,  \
                             ts_zz_xxz,   \
                             ts_zz_xxzz,  \
                             ts_zz_xyy,   \
                             ts_zz_xyyy,  \
                             ts_zz_xyyz,  \
                             ts_zz_xyz,   \
                             ts_zz_xyzz,  \
                             ts_zz_xzz,   \
                             ts_zz_xzzz,  \
                             ts_zz_yyy,   \
                             ts_zz_yyyy,  \
                             ts_zz_yyyz,  \
                             ts_zz_yyz,   \
                             ts_zz_yyzz,  \
                             ts_zz_yzz,   \
                             ts_zz_yzzz,  \
                             ts_zz_zzz,   \
                             ts_zz_zzzz,  \
                             ts_zzz_xxxx, \
                             ts_zzz_xxxy, \
                             ts_zzz_xxxz, \
                             ts_zzz_xxyy, \
                             ts_zzz_xxyz, \
                             ts_zzz_xxzz, \
                             ts_zzz_xyyy, \
                             ts_zzz_xyyz, \
                             ts_zzz_xyzz, \
                             ts_zzz_xzzz, \
                             ts_zzz_yyyy, \
                             ts_zzz_yyyz, \
                             ts_zzz_yyzz, \
                             ts_zzz_yzzz, \
                             ts_zzz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_zzz_xxxx[i] = 2.0 * ts_z_xxxx[i] * fe_0 + ts_zz_xxxx[i] * pa_z[i];

        ts_zzz_xxxy[i] = 2.0 * ts_z_xxxy[i] * fe_0 + ts_zz_xxxy[i] * pa_z[i];

        ts_zzz_xxxz[i] = 2.0 * ts_z_xxxz[i] * fe_0 + ts_zz_xxx[i] * fe_0 + ts_zz_xxxz[i] * pa_z[i];

        ts_zzz_xxyy[i] = 2.0 * ts_z_xxyy[i] * fe_0 + ts_zz_xxyy[i] * pa_z[i];

        ts_zzz_xxyz[i] = 2.0 * ts_z_xxyz[i] * fe_0 + ts_zz_xxy[i] * fe_0 + ts_zz_xxyz[i] * pa_z[i];

        ts_zzz_xxzz[i] = 2.0 * ts_z_xxzz[i] * fe_0 + 2.0 * ts_zz_xxz[i] * fe_0 + ts_zz_xxzz[i] * pa_z[i];

        ts_zzz_xyyy[i] = 2.0 * ts_z_xyyy[i] * fe_0 + ts_zz_xyyy[i] * pa_z[i];

        ts_zzz_xyyz[i] = 2.0 * ts_z_xyyz[i] * fe_0 + ts_zz_xyy[i] * fe_0 + ts_zz_xyyz[i] * pa_z[i];

        ts_zzz_xyzz[i] = 2.0 * ts_z_xyzz[i] * fe_0 + 2.0 * ts_zz_xyz[i] * fe_0 + ts_zz_xyzz[i] * pa_z[i];

        ts_zzz_xzzz[i] = 2.0 * ts_z_xzzz[i] * fe_0 + 3.0 * ts_zz_xzz[i] * fe_0 + ts_zz_xzzz[i] * pa_z[i];

        ts_zzz_yyyy[i] = 2.0 * ts_z_yyyy[i] * fe_0 + ts_zz_yyyy[i] * pa_z[i];

        ts_zzz_yyyz[i] = 2.0 * ts_z_yyyz[i] * fe_0 + ts_zz_yyy[i] * fe_0 + ts_zz_yyyz[i] * pa_z[i];

        ts_zzz_yyzz[i] = 2.0 * ts_z_yyzz[i] * fe_0 + 2.0 * ts_zz_yyz[i] * fe_0 + ts_zz_yyzz[i] * pa_z[i];

        ts_zzz_yzzz[i] = 2.0 * ts_z_yzzz[i] * fe_0 + 3.0 * ts_zz_yzz[i] * fe_0 + ts_zz_yzzz[i] * pa_z[i];

        ts_zzz_zzzz[i] = 2.0 * ts_z_zzzz[i] * fe_0 + 4.0 * ts_zz_zzz[i] * fe_0 + ts_zz_zzzz[i] * pa_z[i];
    }
}

}  // namespace ovlrec
