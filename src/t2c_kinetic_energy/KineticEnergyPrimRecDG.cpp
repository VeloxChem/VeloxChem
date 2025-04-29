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

#include "KineticEnergyPrimRecDG.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_dg(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_dg,
                            const size_t              idx_ovl_sg,
                            const size_t              idx_kin_sg,
                            const size_t              idx_kin_pf,
                            const size_t              idx_kin_pg,
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

    // Set up components of auxiliary buffer : SG

    auto ts_0_xxxx = pbuffer.data(idx_ovl_sg);

    auto ts_0_xxxy = pbuffer.data(idx_ovl_sg + 1);

    auto ts_0_xxxz = pbuffer.data(idx_ovl_sg + 2);

    auto ts_0_xxyy = pbuffer.data(idx_ovl_sg + 3);

    auto ts_0_xxyz = pbuffer.data(idx_ovl_sg + 4);

    auto ts_0_xxzz = pbuffer.data(idx_ovl_sg + 5);

    auto ts_0_xyyy = pbuffer.data(idx_ovl_sg + 6);

    auto ts_0_xyyz = pbuffer.data(idx_ovl_sg + 7);

    auto ts_0_xyzz = pbuffer.data(idx_ovl_sg + 8);

    auto ts_0_xzzz = pbuffer.data(idx_ovl_sg + 9);

    auto ts_0_yyyy = pbuffer.data(idx_ovl_sg + 10);

    auto ts_0_yyyz = pbuffer.data(idx_ovl_sg + 11);

    auto ts_0_yyzz = pbuffer.data(idx_ovl_sg + 12);

    auto ts_0_yzzz = pbuffer.data(idx_ovl_sg + 13);

    auto ts_0_zzzz = pbuffer.data(idx_ovl_sg + 14);

    // Set up components of auxiliary buffer : SG

    auto tk_0_xxxx = pbuffer.data(idx_kin_sg);

    auto tk_0_xxxy = pbuffer.data(idx_kin_sg + 1);

    auto tk_0_xxxz = pbuffer.data(idx_kin_sg + 2);

    auto tk_0_xxyy = pbuffer.data(idx_kin_sg + 3);

    auto tk_0_xxyz = pbuffer.data(idx_kin_sg + 4);

    auto tk_0_xxzz = pbuffer.data(idx_kin_sg + 5);

    auto tk_0_xyyy = pbuffer.data(idx_kin_sg + 6);

    auto tk_0_xyyz = pbuffer.data(idx_kin_sg + 7);

    auto tk_0_xyzz = pbuffer.data(idx_kin_sg + 8);

    auto tk_0_xzzz = pbuffer.data(idx_kin_sg + 9);

    auto tk_0_yyyy = pbuffer.data(idx_kin_sg + 10);

    auto tk_0_yyyz = pbuffer.data(idx_kin_sg + 11);

    auto tk_0_yyzz = pbuffer.data(idx_kin_sg + 12);

    auto tk_0_yzzz = pbuffer.data(idx_kin_sg + 13);

    auto tk_0_zzzz = pbuffer.data(idx_kin_sg + 14);

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

    // Set up components of auxiliary buffer : PG

    auto tk_x_xxxx = pbuffer.data(idx_kin_pg);

    auto tk_x_xxxy = pbuffer.data(idx_kin_pg + 1);

    auto tk_x_xxxz = pbuffer.data(idx_kin_pg + 2);

    auto tk_x_xxyy = pbuffer.data(idx_kin_pg + 3);

    auto tk_x_xxyz = pbuffer.data(idx_kin_pg + 4);

    auto tk_x_xxzz = pbuffer.data(idx_kin_pg + 5);

    auto tk_x_xyyy = pbuffer.data(idx_kin_pg + 6);

    auto tk_x_xyyz = pbuffer.data(idx_kin_pg + 7);

    auto tk_x_xyzz = pbuffer.data(idx_kin_pg + 8);

    auto tk_x_xzzz = pbuffer.data(idx_kin_pg + 9);

    auto tk_x_yyyy = pbuffer.data(idx_kin_pg + 10);

    auto tk_x_yyyz = pbuffer.data(idx_kin_pg + 11);

    auto tk_x_yyzz = pbuffer.data(idx_kin_pg + 12);

    auto tk_x_yzzz = pbuffer.data(idx_kin_pg + 13);

    auto tk_x_zzzz = pbuffer.data(idx_kin_pg + 14);

    auto tk_y_xxxx = pbuffer.data(idx_kin_pg + 15);

    auto tk_y_xxxy = pbuffer.data(idx_kin_pg + 16);

    auto tk_y_xxxz = pbuffer.data(idx_kin_pg + 17);

    auto tk_y_xxyy = pbuffer.data(idx_kin_pg + 18);

    auto tk_y_xxyz = pbuffer.data(idx_kin_pg + 19);

    auto tk_y_xxzz = pbuffer.data(idx_kin_pg + 20);

    auto tk_y_xyyy = pbuffer.data(idx_kin_pg + 21);

    auto tk_y_xyyz = pbuffer.data(idx_kin_pg + 22);

    auto tk_y_xyzz = pbuffer.data(idx_kin_pg + 23);

    auto tk_y_xzzz = pbuffer.data(idx_kin_pg + 24);

    auto tk_y_yyyy = pbuffer.data(idx_kin_pg + 25);

    auto tk_y_yyyz = pbuffer.data(idx_kin_pg + 26);

    auto tk_y_yyzz = pbuffer.data(idx_kin_pg + 27);

    auto tk_y_yzzz = pbuffer.data(idx_kin_pg + 28);

    auto tk_y_zzzz = pbuffer.data(idx_kin_pg + 29);

    auto tk_z_xxxx = pbuffer.data(idx_kin_pg + 30);

    auto tk_z_xxxy = pbuffer.data(idx_kin_pg + 31);

    auto tk_z_xxxz = pbuffer.data(idx_kin_pg + 32);

    auto tk_z_xxyy = pbuffer.data(idx_kin_pg + 33);

    auto tk_z_xxyz = pbuffer.data(idx_kin_pg + 34);

    auto tk_z_xxzz = pbuffer.data(idx_kin_pg + 35);

    auto tk_z_xyyy = pbuffer.data(idx_kin_pg + 36);

    auto tk_z_xyyz = pbuffer.data(idx_kin_pg + 37);

    auto tk_z_xyzz = pbuffer.data(idx_kin_pg + 38);

    auto tk_z_xzzz = pbuffer.data(idx_kin_pg + 39);

    auto tk_z_yyyy = pbuffer.data(idx_kin_pg + 40);

    auto tk_z_yyyz = pbuffer.data(idx_kin_pg + 41);

    auto tk_z_yyzz = pbuffer.data(idx_kin_pg + 42);

    auto tk_z_yzzz = pbuffer.data(idx_kin_pg + 43);

    auto tk_z_zzzz = pbuffer.data(idx_kin_pg + 44);

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

    auto ts_xy_xxxx = pbuffer.data(idx_ovl_dg + 15);

    auto ts_xy_xxxy = pbuffer.data(idx_ovl_dg + 16);

    auto ts_xy_xxxz = pbuffer.data(idx_ovl_dg + 17);

    auto ts_xy_xxyy = pbuffer.data(idx_ovl_dg + 18);

    auto ts_xy_xxyz = pbuffer.data(idx_ovl_dg + 19);

    auto ts_xy_xxzz = pbuffer.data(idx_ovl_dg + 20);

    auto ts_xy_xyyy = pbuffer.data(idx_ovl_dg + 21);

    auto ts_xy_xyyz = pbuffer.data(idx_ovl_dg + 22);

    auto ts_xy_xyzz = pbuffer.data(idx_ovl_dg + 23);

    auto ts_xy_xzzz = pbuffer.data(idx_ovl_dg + 24);

    auto ts_xy_yyyy = pbuffer.data(idx_ovl_dg + 25);

    auto ts_xy_yyyz = pbuffer.data(idx_ovl_dg + 26);

    auto ts_xy_yyzz = pbuffer.data(idx_ovl_dg + 27);

    auto ts_xy_yzzz = pbuffer.data(idx_ovl_dg + 28);

    auto ts_xy_zzzz = pbuffer.data(idx_ovl_dg + 29);

    auto ts_xz_xxxx = pbuffer.data(idx_ovl_dg + 30);

    auto ts_xz_xxxy = pbuffer.data(idx_ovl_dg + 31);

    auto ts_xz_xxxz = pbuffer.data(idx_ovl_dg + 32);

    auto ts_xz_xxyy = pbuffer.data(idx_ovl_dg + 33);

    auto ts_xz_xxyz = pbuffer.data(idx_ovl_dg + 34);

    auto ts_xz_xxzz = pbuffer.data(idx_ovl_dg + 35);

    auto ts_xz_xyyy = pbuffer.data(idx_ovl_dg + 36);

    auto ts_xz_xyyz = pbuffer.data(idx_ovl_dg + 37);

    auto ts_xz_xyzz = pbuffer.data(idx_ovl_dg + 38);

    auto ts_xz_xzzz = pbuffer.data(idx_ovl_dg + 39);

    auto ts_xz_yyyy = pbuffer.data(idx_ovl_dg + 40);

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

    auto ts_yz_xxxx = pbuffer.data(idx_ovl_dg + 60);

    auto ts_yz_xxxy = pbuffer.data(idx_ovl_dg + 61);

    auto ts_yz_xxxz = pbuffer.data(idx_ovl_dg + 62);

    auto ts_yz_xxyy = pbuffer.data(idx_ovl_dg + 63);

    auto ts_yz_xxyz = pbuffer.data(idx_ovl_dg + 64);

    auto ts_yz_xxzz = pbuffer.data(idx_ovl_dg + 65);

    auto ts_yz_xyyy = pbuffer.data(idx_ovl_dg + 66);

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

    // Set up 0-15 components of targeted buffer : DG

    auto tk_xx_xxxx = pbuffer.data(idx_kin_dg);

    auto tk_xx_xxxy = pbuffer.data(idx_kin_dg + 1);

    auto tk_xx_xxxz = pbuffer.data(idx_kin_dg + 2);

    auto tk_xx_xxyy = pbuffer.data(idx_kin_dg + 3);

    auto tk_xx_xxyz = pbuffer.data(idx_kin_dg + 4);

    auto tk_xx_xxzz = pbuffer.data(idx_kin_dg + 5);

    auto tk_xx_xyyy = pbuffer.data(idx_kin_dg + 6);

    auto tk_xx_xyyz = pbuffer.data(idx_kin_dg + 7);

    auto tk_xx_xyzz = pbuffer.data(idx_kin_dg + 8);

    auto tk_xx_xzzz = pbuffer.data(idx_kin_dg + 9);

    auto tk_xx_yyyy = pbuffer.data(idx_kin_dg + 10);

    auto tk_xx_yyyz = pbuffer.data(idx_kin_dg + 11);

    auto tk_xx_yyzz = pbuffer.data(idx_kin_dg + 12);

    auto tk_xx_yzzz = pbuffer.data(idx_kin_dg + 13);

    auto tk_xx_zzzz = pbuffer.data(idx_kin_dg + 14);

#pragma omp simd aligned(pa_x,           \
                             tk_0_xxxx,  \
                             tk_0_xxxy,  \
                             tk_0_xxxz,  \
                             tk_0_xxyy,  \
                             tk_0_xxyz,  \
                             tk_0_xxzz,  \
                             tk_0_xyyy,  \
                             tk_0_xyyz,  \
                             tk_0_xyzz,  \
                             tk_0_xzzz,  \
                             tk_0_yyyy,  \
                             tk_0_yyyz,  \
                             tk_0_yyzz,  \
                             tk_0_yzzz,  \
                             tk_0_zzzz,  \
                             tk_x_xxx,   \
                             tk_x_xxxx,  \
                             tk_x_xxxy,  \
                             tk_x_xxxz,  \
                             tk_x_xxy,   \
                             tk_x_xxyy,  \
                             tk_x_xxyz,  \
                             tk_x_xxz,   \
                             tk_x_xxzz,  \
                             tk_x_xyy,   \
                             tk_x_xyyy,  \
                             tk_x_xyyz,  \
                             tk_x_xyz,   \
                             tk_x_xyzz,  \
                             tk_x_xzz,   \
                             tk_x_xzzz,  \
                             tk_x_yyy,   \
                             tk_x_yyyy,  \
                             tk_x_yyyz,  \
                             tk_x_yyz,   \
                             tk_x_yyzz,  \
                             tk_x_yzz,   \
                             tk_x_yzzz,  \
                             tk_x_zzz,   \
                             tk_x_zzzz,  \
                             tk_xx_xxxx, \
                             tk_xx_xxxy, \
                             tk_xx_xxxz, \
                             tk_xx_xxyy, \
                             tk_xx_xxyz, \
                             tk_xx_xxzz, \
                             tk_xx_xyyy, \
                             tk_xx_xyyz, \
                             tk_xx_xyzz, \
                             tk_xx_xzzz, \
                             tk_xx_yyyy, \
                             tk_xx_yyyz, \
                             tk_xx_yyzz, \
                             tk_xx_yzzz, \
                             tk_xx_zzzz, \
                             ts_0_xxxx,  \
                             ts_0_xxxy,  \
                             ts_0_xxxz,  \
                             ts_0_xxyy,  \
                             ts_0_xxyz,  \
                             ts_0_xxzz,  \
                             ts_0_xyyy,  \
                             ts_0_xyyz,  \
                             ts_0_xyzz,  \
                             ts_0_xzzz,  \
                             ts_0_yyyy,  \
                             ts_0_yyyz,  \
                             ts_0_yyzz,  \
                             ts_0_yzzz,  \
                             ts_0_zzzz,  \
                             ts_xx_xxxx, \
                             ts_xx_xxxy, \
                             ts_xx_xxxz, \
                             ts_xx_xxyy, \
                             ts_xx_xxyz, \
                             ts_xx_xxzz, \
                             ts_xx_xyyy, \
                             ts_xx_xyyz, \
                             ts_xx_xyzz, \
                             ts_xx_xzzz, \
                             ts_xx_yyyy, \
                             ts_xx_yyyz, \
                             ts_xx_yyzz, \
                             ts_xx_yzzz, \
                             ts_xx_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xx_xxxx[i] =
            -2.0 * ts_0_xxxx[i] * fbe_0 * fz_0 + tk_0_xxxx[i] * fe_0 + 4.0 * tk_x_xxx[i] * fe_0 + tk_x_xxxx[i] * pa_x[i] + 2.0 * ts_xx_xxxx[i] * fz_0;

        tk_xx_xxxy[i] =
            -2.0 * ts_0_xxxy[i] * fbe_0 * fz_0 + tk_0_xxxy[i] * fe_0 + 3.0 * tk_x_xxy[i] * fe_0 + tk_x_xxxy[i] * pa_x[i] + 2.0 * ts_xx_xxxy[i] * fz_0;

        tk_xx_xxxz[i] =
            -2.0 * ts_0_xxxz[i] * fbe_0 * fz_0 + tk_0_xxxz[i] * fe_0 + 3.0 * tk_x_xxz[i] * fe_0 + tk_x_xxxz[i] * pa_x[i] + 2.0 * ts_xx_xxxz[i] * fz_0;

        tk_xx_xxyy[i] =
            -2.0 * ts_0_xxyy[i] * fbe_0 * fz_0 + tk_0_xxyy[i] * fe_0 + 2.0 * tk_x_xyy[i] * fe_0 + tk_x_xxyy[i] * pa_x[i] + 2.0 * ts_xx_xxyy[i] * fz_0;

        tk_xx_xxyz[i] =
            -2.0 * ts_0_xxyz[i] * fbe_0 * fz_0 + tk_0_xxyz[i] * fe_0 + 2.0 * tk_x_xyz[i] * fe_0 + tk_x_xxyz[i] * pa_x[i] + 2.0 * ts_xx_xxyz[i] * fz_0;

        tk_xx_xxzz[i] =
            -2.0 * ts_0_xxzz[i] * fbe_0 * fz_0 + tk_0_xxzz[i] * fe_0 + 2.0 * tk_x_xzz[i] * fe_0 + tk_x_xxzz[i] * pa_x[i] + 2.0 * ts_xx_xxzz[i] * fz_0;

        tk_xx_xyyy[i] =
            -2.0 * ts_0_xyyy[i] * fbe_0 * fz_0 + tk_0_xyyy[i] * fe_0 + tk_x_yyy[i] * fe_0 + tk_x_xyyy[i] * pa_x[i] + 2.0 * ts_xx_xyyy[i] * fz_0;

        tk_xx_xyyz[i] =
            -2.0 * ts_0_xyyz[i] * fbe_0 * fz_0 + tk_0_xyyz[i] * fe_0 + tk_x_yyz[i] * fe_0 + tk_x_xyyz[i] * pa_x[i] + 2.0 * ts_xx_xyyz[i] * fz_0;

        tk_xx_xyzz[i] =
            -2.0 * ts_0_xyzz[i] * fbe_0 * fz_0 + tk_0_xyzz[i] * fe_0 + tk_x_yzz[i] * fe_0 + tk_x_xyzz[i] * pa_x[i] + 2.0 * ts_xx_xyzz[i] * fz_0;

        tk_xx_xzzz[i] =
            -2.0 * ts_0_xzzz[i] * fbe_0 * fz_0 + tk_0_xzzz[i] * fe_0 + tk_x_zzz[i] * fe_0 + tk_x_xzzz[i] * pa_x[i] + 2.0 * ts_xx_xzzz[i] * fz_0;

        tk_xx_yyyy[i] = -2.0 * ts_0_yyyy[i] * fbe_0 * fz_0 + tk_0_yyyy[i] * fe_0 + tk_x_yyyy[i] * pa_x[i] + 2.0 * ts_xx_yyyy[i] * fz_0;

        tk_xx_yyyz[i] = -2.0 * ts_0_yyyz[i] * fbe_0 * fz_0 + tk_0_yyyz[i] * fe_0 + tk_x_yyyz[i] * pa_x[i] + 2.0 * ts_xx_yyyz[i] * fz_0;

        tk_xx_yyzz[i] = -2.0 * ts_0_yyzz[i] * fbe_0 * fz_0 + tk_0_yyzz[i] * fe_0 + tk_x_yyzz[i] * pa_x[i] + 2.0 * ts_xx_yyzz[i] * fz_0;

        tk_xx_yzzz[i] = -2.0 * ts_0_yzzz[i] * fbe_0 * fz_0 + tk_0_yzzz[i] * fe_0 + tk_x_yzzz[i] * pa_x[i] + 2.0 * ts_xx_yzzz[i] * fz_0;

        tk_xx_zzzz[i] = -2.0 * ts_0_zzzz[i] * fbe_0 * fz_0 + tk_0_zzzz[i] * fe_0 + tk_x_zzzz[i] * pa_x[i] + 2.0 * ts_xx_zzzz[i] * fz_0;
    }

    // Set up 15-30 components of targeted buffer : DG

    auto tk_xy_xxxx = pbuffer.data(idx_kin_dg + 15);

    auto tk_xy_xxxy = pbuffer.data(idx_kin_dg + 16);

    auto tk_xy_xxxz = pbuffer.data(idx_kin_dg + 17);

    auto tk_xy_xxyy = pbuffer.data(idx_kin_dg + 18);

    auto tk_xy_xxyz = pbuffer.data(idx_kin_dg + 19);

    auto tk_xy_xxzz = pbuffer.data(idx_kin_dg + 20);

    auto tk_xy_xyyy = pbuffer.data(idx_kin_dg + 21);

    auto tk_xy_xyyz = pbuffer.data(idx_kin_dg + 22);

    auto tk_xy_xyzz = pbuffer.data(idx_kin_dg + 23);

    auto tk_xy_xzzz = pbuffer.data(idx_kin_dg + 24);

    auto tk_xy_yyyy = pbuffer.data(idx_kin_dg + 25);

    auto tk_xy_yyyz = pbuffer.data(idx_kin_dg + 26);

    auto tk_xy_yyzz = pbuffer.data(idx_kin_dg + 27);

    auto tk_xy_yzzz = pbuffer.data(idx_kin_dg + 28);

    auto tk_xy_zzzz = pbuffer.data(idx_kin_dg + 29);

#pragma omp simd aligned(pa_x,           \
                             pa_y,       \
                             tk_x_xxxx,  \
                             tk_x_xxxz,  \
                             tk_x_xxzz,  \
                             tk_x_xzzz,  \
                             tk_xy_xxxx, \
                             tk_xy_xxxy, \
                             tk_xy_xxxz, \
                             tk_xy_xxyy, \
                             tk_xy_xxyz, \
                             tk_xy_xxzz, \
                             tk_xy_xyyy, \
                             tk_xy_xyyz, \
                             tk_xy_xyzz, \
                             tk_xy_xzzz, \
                             tk_xy_yyyy, \
                             tk_xy_yyyz, \
                             tk_xy_yyzz, \
                             tk_xy_yzzz, \
                             tk_xy_zzzz, \
                             tk_y_xxxy,  \
                             tk_y_xxy,   \
                             tk_y_xxyy,  \
                             tk_y_xxyz,  \
                             tk_y_xyy,   \
                             tk_y_xyyy,  \
                             tk_y_xyyz,  \
                             tk_y_xyz,   \
                             tk_y_xyzz,  \
                             tk_y_yyy,   \
                             tk_y_yyyy,  \
                             tk_y_yyyz,  \
                             tk_y_yyz,   \
                             tk_y_yyzz,  \
                             tk_y_yzz,   \
                             tk_y_yzzz,  \
                             tk_y_zzzz,  \
                             ts_xy_xxxx, \
                             ts_xy_xxxy, \
                             ts_xy_xxxz, \
                             ts_xy_xxyy, \
                             ts_xy_xxyz, \
                             ts_xy_xxzz, \
                             ts_xy_xyyy, \
                             ts_xy_xyyz, \
                             ts_xy_xyzz, \
                             ts_xy_xzzz, \
                             ts_xy_yyyy, \
                             ts_xy_yyyz, \
                             ts_xy_yyzz, \
                             ts_xy_yzzz, \
                             ts_xy_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xy_xxxx[i] = tk_x_xxxx[i] * pa_y[i] + 2.0 * ts_xy_xxxx[i] * fz_0;

        tk_xy_xxxy[i] = 3.0 * tk_y_xxy[i] * fe_0 + tk_y_xxxy[i] * pa_x[i] + 2.0 * ts_xy_xxxy[i] * fz_0;

        tk_xy_xxxz[i] = tk_x_xxxz[i] * pa_y[i] + 2.0 * ts_xy_xxxz[i] * fz_0;

        tk_xy_xxyy[i] = 2.0 * tk_y_xyy[i] * fe_0 + tk_y_xxyy[i] * pa_x[i] + 2.0 * ts_xy_xxyy[i] * fz_0;

        tk_xy_xxyz[i] = 2.0 * tk_y_xyz[i] * fe_0 + tk_y_xxyz[i] * pa_x[i] + 2.0 * ts_xy_xxyz[i] * fz_0;

        tk_xy_xxzz[i] = tk_x_xxzz[i] * pa_y[i] + 2.0 * ts_xy_xxzz[i] * fz_0;

        tk_xy_xyyy[i] = tk_y_yyy[i] * fe_0 + tk_y_xyyy[i] * pa_x[i] + 2.0 * ts_xy_xyyy[i] * fz_0;

        tk_xy_xyyz[i] = tk_y_yyz[i] * fe_0 + tk_y_xyyz[i] * pa_x[i] + 2.0 * ts_xy_xyyz[i] * fz_0;

        tk_xy_xyzz[i] = tk_y_yzz[i] * fe_0 + tk_y_xyzz[i] * pa_x[i] + 2.0 * ts_xy_xyzz[i] * fz_0;

        tk_xy_xzzz[i] = tk_x_xzzz[i] * pa_y[i] + 2.0 * ts_xy_xzzz[i] * fz_0;

        tk_xy_yyyy[i] = tk_y_yyyy[i] * pa_x[i] + 2.0 * ts_xy_yyyy[i] * fz_0;

        tk_xy_yyyz[i] = tk_y_yyyz[i] * pa_x[i] + 2.0 * ts_xy_yyyz[i] * fz_0;

        tk_xy_yyzz[i] = tk_y_yyzz[i] * pa_x[i] + 2.0 * ts_xy_yyzz[i] * fz_0;

        tk_xy_yzzz[i] = tk_y_yzzz[i] * pa_x[i] + 2.0 * ts_xy_yzzz[i] * fz_0;

        tk_xy_zzzz[i] = tk_y_zzzz[i] * pa_x[i] + 2.0 * ts_xy_zzzz[i] * fz_0;
    }

    // Set up 30-45 components of targeted buffer : DG

    auto tk_xz_xxxx = pbuffer.data(idx_kin_dg + 30);

    auto tk_xz_xxxy = pbuffer.data(idx_kin_dg + 31);

    auto tk_xz_xxxz = pbuffer.data(idx_kin_dg + 32);

    auto tk_xz_xxyy = pbuffer.data(idx_kin_dg + 33);

    auto tk_xz_xxyz = pbuffer.data(idx_kin_dg + 34);

    auto tk_xz_xxzz = pbuffer.data(idx_kin_dg + 35);

    auto tk_xz_xyyy = pbuffer.data(idx_kin_dg + 36);

    auto tk_xz_xyyz = pbuffer.data(idx_kin_dg + 37);

    auto tk_xz_xyzz = pbuffer.data(idx_kin_dg + 38);

    auto tk_xz_xzzz = pbuffer.data(idx_kin_dg + 39);

    auto tk_xz_yyyy = pbuffer.data(idx_kin_dg + 40);

    auto tk_xz_yyyz = pbuffer.data(idx_kin_dg + 41);

    auto tk_xz_yyzz = pbuffer.data(idx_kin_dg + 42);

    auto tk_xz_yzzz = pbuffer.data(idx_kin_dg + 43);

    auto tk_xz_zzzz = pbuffer.data(idx_kin_dg + 44);

#pragma omp simd aligned(pa_x,           \
                             pa_z,       \
                             tk_x_xxxx,  \
                             tk_x_xxxy,  \
                             tk_x_xxyy,  \
                             tk_x_xyyy,  \
                             tk_xz_xxxx, \
                             tk_xz_xxxy, \
                             tk_xz_xxxz, \
                             tk_xz_xxyy, \
                             tk_xz_xxyz, \
                             tk_xz_xxzz, \
                             tk_xz_xyyy, \
                             tk_xz_xyyz, \
                             tk_xz_xyzz, \
                             tk_xz_xzzz, \
                             tk_xz_yyyy, \
                             tk_xz_yyyz, \
                             tk_xz_yyzz, \
                             tk_xz_yzzz, \
                             tk_xz_zzzz, \
                             tk_z_xxxz,  \
                             tk_z_xxyz,  \
                             tk_z_xxz,   \
                             tk_z_xxzz,  \
                             tk_z_xyyz,  \
                             tk_z_xyz,   \
                             tk_z_xyzz,  \
                             tk_z_xzz,   \
                             tk_z_xzzz,  \
                             tk_z_yyyy,  \
                             tk_z_yyyz,  \
                             tk_z_yyz,   \
                             tk_z_yyzz,  \
                             tk_z_yzz,   \
                             tk_z_yzzz,  \
                             tk_z_zzz,   \
                             tk_z_zzzz,  \
                             ts_xz_xxxx, \
                             ts_xz_xxxy, \
                             ts_xz_xxxz, \
                             ts_xz_xxyy, \
                             ts_xz_xxyz, \
                             ts_xz_xxzz, \
                             ts_xz_xyyy, \
                             ts_xz_xyyz, \
                             ts_xz_xyzz, \
                             ts_xz_xzzz, \
                             ts_xz_yyyy, \
                             ts_xz_yyyz, \
                             ts_xz_yyzz, \
                             ts_xz_yzzz, \
                             ts_xz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xz_xxxx[i] = tk_x_xxxx[i] * pa_z[i] + 2.0 * ts_xz_xxxx[i] * fz_0;

        tk_xz_xxxy[i] = tk_x_xxxy[i] * pa_z[i] + 2.0 * ts_xz_xxxy[i] * fz_0;

        tk_xz_xxxz[i] = 3.0 * tk_z_xxz[i] * fe_0 + tk_z_xxxz[i] * pa_x[i] + 2.0 * ts_xz_xxxz[i] * fz_0;

        tk_xz_xxyy[i] = tk_x_xxyy[i] * pa_z[i] + 2.0 * ts_xz_xxyy[i] * fz_0;

        tk_xz_xxyz[i] = 2.0 * tk_z_xyz[i] * fe_0 + tk_z_xxyz[i] * pa_x[i] + 2.0 * ts_xz_xxyz[i] * fz_0;

        tk_xz_xxzz[i] = 2.0 * tk_z_xzz[i] * fe_0 + tk_z_xxzz[i] * pa_x[i] + 2.0 * ts_xz_xxzz[i] * fz_0;

        tk_xz_xyyy[i] = tk_x_xyyy[i] * pa_z[i] + 2.0 * ts_xz_xyyy[i] * fz_0;

        tk_xz_xyyz[i] = tk_z_yyz[i] * fe_0 + tk_z_xyyz[i] * pa_x[i] + 2.0 * ts_xz_xyyz[i] * fz_0;

        tk_xz_xyzz[i] = tk_z_yzz[i] * fe_0 + tk_z_xyzz[i] * pa_x[i] + 2.0 * ts_xz_xyzz[i] * fz_0;

        tk_xz_xzzz[i] = tk_z_zzz[i] * fe_0 + tk_z_xzzz[i] * pa_x[i] + 2.0 * ts_xz_xzzz[i] * fz_0;

        tk_xz_yyyy[i] = tk_z_yyyy[i] * pa_x[i] + 2.0 * ts_xz_yyyy[i] * fz_0;

        tk_xz_yyyz[i] = tk_z_yyyz[i] * pa_x[i] + 2.0 * ts_xz_yyyz[i] * fz_0;

        tk_xz_yyzz[i] = tk_z_yyzz[i] * pa_x[i] + 2.0 * ts_xz_yyzz[i] * fz_0;

        tk_xz_yzzz[i] = tk_z_yzzz[i] * pa_x[i] + 2.0 * ts_xz_yzzz[i] * fz_0;

        tk_xz_zzzz[i] = tk_z_zzzz[i] * pa_x[i] + 2.0 * ts_xz_zzzz[i] * fz_0;
    }

    // Set up 45-60 components of targeted buffer : DG

    auto tk_yy_xxxx = pbuffer.data(idx_kin_dg + 45);

    auto tk_yy_xxxy = pbuffer.data(idx_kin_dg + 46);

    auto tk_yy_xxxz = pbuffer.data(idx_kin_dg + 47);

    auto tk_yy_xxyy = pbuffer.data(idx_kin_dg + 48);

    auto tk_yy_xxyz = pbuffer.data(idx_kin_dg + 49);

    auto tk_yy_xxzz = pbuffer.data(idx_kin_dg + 50);

    auto tk_yy_xyyy = pbuffer.data(idx_kin_dg + 51);

    auto tk_yy_xyyz = pbuffer.data(idx_kin_dg + 52);

    auto tk_yy_xyzz = pbuffer.data(idx_kin_dg + 53);

    auto tk_yy_xzzz = pbuffer.data(idx_kin_dg + 54);

    auto tk_yy_yyyy = pbuffer.data(idx_kin_dg + 55);

    auto tk_yy_yyyz = pbuffer.data(idx_kin_dg + 56);

    auto tk_yy_yyzz = pbuffer.data(idx_kin_dg + 57);

    auto tk_yy_yzzz = pbuffer.data(idx_kin_dg + 58);

    auto tk_yy_zzzz = pbuffer.data(idx_kin_dg + 59);

#pragma omp simd aligned(pa_y,           \
                             tk_0_xxxx,  \
                             tk_0_xxxy,  \
                             tk_0_xxxz,  \
                             tk_0_xxyy,  \
                             tk_0_xxyz,  \
                             tk_0_xxzz,  \
                             tk_0_xyyy,  \
                             tk_0_xyyz,  \
                             tk_0_xyzz,  \
                             tk_0_xzzz,  \
                             tk_0_yyyy,  \
                             tk_0_yyyz,  \
                             tk_0_yyzz,  \
                             tk_0_yzzz,  \
                             tk_0_zzzz,  \
                             tk_y_xxx,   \
                             tk_y_xxxx,  \
                             tk_y_xxxy,  \
                             tk_y_xxxz,  \
                             tk_y_xxy,   \
                             tk_y_xxyy,  \
                             tk_y_xxyz,  \
                             tk_y_xxz,   \
                             tk_y_xxzz,  \
                             tk_y_xyy,   \
                             tk_y_xyyy,  \
                             tk_y_xyyz,  \
                             tk_y_xyz,   \
                             tk_y_xyzz,  \
                             tk_y_xzz,   \
                             tk_y_xzzz,  \
                             tk_y_yyy,   \
                             tk_y_yyyy,  \
                             tk_y_yyyz,  \
                             tk_y_yyz,   \
                             tk_y_yyzz,  \
                             tk_y_yzz,   \
                             tk_y_yzzz,  \
                             tk_y_zzz,   \
                             tk_y_zzzz,  \
                             tk_yy_xxxx, \
                             tk_yy_xxxy, \
                             tk_yy_xxxz, \
                             tk_yy_xxyy, \
                             tk_yy_xxyz, \
                             tk_yy_xxzz, \
                             tk_yy_xyyy, \
                             tk_yy_xyyz, \
                             tk_yy_xyzz, \
                             tk_yy_xzzz, \
                             tk_yy_yyyy, \
                             tk_yy_yyyz, \
                             tk_yy_yyzz, \
                             tk_yy_yzzz, \
                             tk_yy_zzzz, \
                             ts_0_xxxx,  \
                             ts_0_xxxy,  \
                             ts_0_xxxz,  \
                             ts_0_xxyy,  \
                             ts_0_xxyz,  \
                             ts_0_xxzz,  \
                             ts_0_xyyy,  \
                             ts_0_xyyz,  \
                             ts_0_xyzz,  \
                             ts_0_xzzz,  \
                             ts_0_yyyy,  \
                             ts_0_yyyz,  \
                             ts_0_yyzz,  \
                             ts_0_yzzz,  \
                             ts_0_zzzz,  \
                             ts_yy_xxxx, \
                             ts_yy_xxxy, \
                             ts_yy_xxxz, \
                             ts_yy_xxyy, \
                             ts_yy_xxyz, \
                             ts_yy_xxzz, \
                             ts_yy_xyyy, \
                             ts_yy_xyyz, \
                             ts_yy_xyzz, \
                             ts_yy_xzzz, \
                             ts_yy_yyyy, \
                             ts_yy_yyyz, \
                             ts_yy_yyzz, \
                             ts_yy_yzzz, \
                             ts_yy_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yy_xxxx[i] = -2.0 * ts_0_xxxx[i] * fbe_0 * fz_0 + tk_0_xxxx[i] * fe_0 + tk_y_xxxx[i] * pa_y[i] + 2.0 * ts_yy_xxxx[i] * fz_0;

        tk_yy_xxxy[i] =
            -2.0 * ts_0_xxxy[i] * fbe_0 * fz_0 + tk_0_xxxy[i] * fe_0 + tk_y_xxx[i] * fe_0 + tk_y_xxxy[i] * pa_y[i] + 2.0 * ts_yy_xxxy[i] * fz_0;

        tk_yy_xxxz[i] = -2.0 * ts_0_xxxz[i] * fbe_0 * fz_0 + tk_0_xxxz[i] * fe_0 + tk_y_xxxz[i] * pa_y[i] + 2.0 * ts_yy_xxxz[i] * fz_0;

        tk_yy_xxyy[i] =
            -2.0 * ts_0_xxyy[i] * fbe_0 * fz_0 + tk_0_xxyy[i] * fe_0 + 2.0 * tk_y_xxy[i] * fe_0 + tk_y_xxyy[i] * pa_y[i] + 2.0 * ts_yy_xxyy[i] * fz_0;

        tk_yy_xxyz[i] =
            -2.0 * ts_0_xxyz[i] * fbe_0 * fz_0 + tk_0_xxyz[i] * fe_0 + tk_y_xxz[i] * fe_0 + tk_y_xxyz[i] * pa_y[i] + 2.0 * ts_yy_xxyz[i] * fz_0;

        tk_yy_xxzz[i] = -2.0 * ts_0_xxzz[i] * fbe_0 * fz_0 + tk_0_xxzz[i] * fe_0 + tk_y_xxzz[i] * pa_y[i] + 2.0 * ts_yy_xxzz[i] * fz_0;

        tk_yy_xyyy[i] =
            -2.0 * ts_0_xyyy[i] * fbe_0 * fz_0 + tk_0_xyyy[i] * fe_0 + 3.0 * tk_y_xyy[i] * fe_0 + tk_y_xyyy[i] * pa_y[i] + 2.0 * ts_yy_xyyy[i] * fz_0;

        tk_yy_xyyz[i] =
            -2.0 * ts_0_xyyz[i] * fbe_0 * fz_0 + tk_0_xyyz[i] * fe_0 + 2.0 * tk_y_xyz[i] * fe_0 + tk_y_xyyz[i] * pa_y[i] + 2.0 * ts_yy_xyyz[i] * fz_0;

        tk_yy_xyzz[i] =
            -2.0 * ts_0_xyzz[i] * fbe_0 * fz_0 + tk_0_xyzz[i] * fe_0 + tk_y_xzz[i] * fe_0 + tk_y_xyzz[i] * pa_y[i] + 2.0 * ts_yy_xyzz[i] * fz_0;

        tk_yy_xzzz[i] = -2.0 * ts_0_xzzz[i] * fbe_0 * fz_0 + tk_0_xzzz[i] * fe_0 + tk_y_xzzz[i] * pa_y[i] + 2.0 * ts_yy_xzzz[i] * fz_0;

        tk_yy_yyyy[i] =
            -2.0 * ts_0_yyyy[i] * fbe_0 * fz_0 + tk_0_yyyy[i] * fe_0 + 4.0 * tk_y_yyy[i] * fe_0 + tk_y_yyyy[i] * pa_y[i] + 2.0 * ts_yy_yyyy[i] * fz_0;

        tk_yy_yyyz[i] =
            -2.0 * ts_0_yyyz[i] * fbe_0 * fz_0 + tk_0_yyyz[i] * fe_0 + 3.0 * tk_y_yyz[i] * fe_0 + tk_y_yyyz[i] * pa_y[i] + 2.0 * ts_yy_yyyz[i] * fz_0;

        tk_yy_yyzz[i] =
            -2.0 * ts_0_yyzz[i] * fbe_0 * fz_0 + tk_0_yyzz[i] * fe_0 + 2.0 * tk_y_yzz[i] * fe_0 + tk_y_yyzz[i] * pa_y[i] + 2.0 * ts_yy_yyzz[i] * fz_0;

        tk_yy_yzzz[i] =
            -2.0 * ts_0_yzzz[i] * fbe_0 * fz_0 + tk_0_yzzz[i] * fe_0 + tk_y_zzz[i] * fe_0 + tk_y_yzzz[i] * pa_y[i] + 2.0 * ts_yy_yzzz[i] * fz_0;

        tk_yy_zzzz[i] = -2.0 * ts_0_zzzz[i] * fbe_0 * fz_0 + tk_0_zzzz[i] * fe_0 + tk_y_zzzz[i] * pa_y[i] + 2.0 * ts_yy_zzzz[i] * fz_0;
    }

    // Set up 60-75 components of targeted buffer : DG

    auto tk_yz_xxxx = pbuffer.data(idx_kin_dg + 60);

    auto tk_yz_xxxy = pbuffer.data(idx_kin_dg + 61);

    auto tk_yz_xxxz = pbuffer.data(idx_kin_dg + 62);

    auto tk_yz_xxyy = pbuffer.data(idx_kin_dg + 63);

    auto tk_yz_xxyz = pbuffer.data(idx_kin_dg + 64);

    auto tk_yz_xxzz = pbuffer.data(idx_kin_dg + 65);

    auto tk_yz_xyyy = pbuffer.data(idx_kin_dg + 66);

    auto tk_yz_xyyz = pbuffer.data(idx_kin_dg + 67);

    auto tk_yz_xyzz = pbuffer.data(idx_kin_dg + 68);

    auto tk_yz_xzzz = pbuffer.data(idx_kin_dg + 69);

    auto tk_yz_yyyy = pbuffer.data(idx_kin_dg + 70);

    auto tk_yz_yyyz = pbuffer.data(idx_kin_dg + 71);

    auto tk_yz_yyzz = pbuffer.data(idx_kin_dg + 72);

    auto tk_yz_yzzz = pbuffer.data(idx_kin_dg + 73);

    auto tk_yz_zzzz = pbuffer.data(idx_kin_dg + 74);

#pragma omp simd aligned(pa_y,           \
                             pa_z,       \
                             tk_y_xxxy,  \
                             tk_y_xxyy,  \
                             tk_y_xyyy,  \
                             tk_y_yyyy,  \
                             tk_yz_xxxx, \
                             tk_yz_xxxy, \
                             tk_yz_xxxz, \
                             tk_yz_xxyy, \
                             tk_yz_xxyz, \
                             tk_yz_xxzz, \
                             tk_yz_xyyy, \
                             tk_yz_xyyz, \
                             tk_yz_xyzz, \
                             tk_yz_xzzz, \
                             tk_yz_yyyy, \
                             tk_yz_yyyz, \
                             tk_yz_yyzz, \
                             tk_yz_yzzz, \
                             tk_yz_zzzz, \
                             tk_z_xxxx,  \
                             tk_z_xxxz,  \
                             tk_z_xxyz,  \
                             tk_z_xxz,   \
                             tk_z_xxzz,  \
                             tk_z_xyyz,  \
                             tk_z_xyz,   \
                             tk_z_xyzz,  \
                             tk_z_xzz,   \
                             tk_z_xzzz,  \
                             tk_z_yyyz,  \
                             tk_z_yyz,   \
                             tk_z_yyzz,  \
                             tk_z_yzz,   \
                             tk_z_yzzz,  \
                             tk_z_zzz,   \
                             tk_z_zzzz,  \
                             ts_yz_xxxx, \
                             ts_yz_xxxy, \
                             ts_yz_xxxz, \
                             ts_yz_xxyy, \
                             ts_yz_xxyz, \
                             ts_yz_xxzz, \
                             ts_yz_xyyy, \
                             ts_yz_xyyz, \
                             ts_yz_xyzz, \
                             ts_yz_xzzz, \
                             ts_yz_yyyy, \
                             ts_yz_yyyz, \
                             ts_yz_yyzz, \
                             ts_yz_yzzz, \
                             ts_yz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yz_xxxx[i] = tk_z_xxxx[i] * pa_y[i] + 2.0 * ts_yz_xxxx[i] * fz_0;

        tk_yz_xxxy[i] = tk_y_xxxy[i] * pa_z[i] + 2.0 * ts_yz_xxxy[i] * fz_0;

        tk_yz_xxxz[i] = tk_z_xxxz[i] * pa_y[i] + 2.0 * ts_yz_xxxz[i] * fz_0;

        tk_yz_xxyy[i] = tk_y_xxyy[i] * pa_z[i] + 2.0 * ts_yz_xxyy[i] * fz_0;

        tk_yz_xxyz[i] = tk_z_xxz[i] * fe_0 + tk_z_xxyz[i] * pa_y[i] + 2.0 * ts_yz_xxyz[i] * fz_0;

        tk_yz_xxzz[i] = tk_z_xxzz[i] * pa_y[i] + 2.0 * ts_yz_xxzz[i] * fz_0;

        tk_yz_xyyy[i] = tk_y_xyyy[i] * pa_z[i] + 2.0 * ts_yz_xyyy[i] * fz_0;

        tk_yz_xyyz[i] = 2.0 * tk_z_xyz[i] * fe_0 + tk_z_xyyz[i] * pa_y[i] + 2.0 * ts_yz_xyyz[i] * fz_0;

        tk_yz_xyzz[i] = tk_z_xzz[i] * fe_0 + tk_z_xyzz[i] * pa_y[i] + 2.0 * ts_yz_xyzz[i] * fz_0;

        tk_yz_xzzz[i] = tk_z_xzzz[i] * pa_y[i] + 2.0 * ts_yz_xzzz[i] * fz_0;

        tk_yz_yyyy[i] = tk_y_yyyy[i] * pa_z[i] + 2.0 * ts_yz_yyyy[i] * fz_0;

        tk_yz_yyyz[i] = 3.0 * tk_z_yyz[i] * fe_0 + tk_z_yyyz[i] * pa_y[i] + 2.0 * ts_yz_yyyz[i] * fz_0;

        tk_yz_yyzz[i] = 2.0 * tk_z_yzz[i] * fe_0 + tk_z_yyzz[i] * pa_y[i] + 2.0 * ts_yz_yyzz[i] * fz_0;

        tk_yz_yzzz[i] = tk_z_zzz[i] * fe_0 + tk_z_yzzz[i] * pa_y[i] + 2.0 * ts_yz_yzzz[i] * fz_0;

        tk_yz_zzzz[i] = tk_z_zzzz[i] * pa_y[i] + 2.0 * ts_yz_zzzz[i] * fz_0;
    }

    // Set up 75-90 components of targeted buffer : DG

    auto tk_zz_xxxx = pbuffer.data(idx_kin_dg + 75);

    auto tk_zz_xxxy = pbuffer.data(idx_kin_dg + 76);

    auto tk_zz_xxxz = pbuffer.data(idx_kin_dg + 77);

    auto tk_zz_xxyy = pbuffer.data(idx_kin_dg + 78);

    auto tk_zz_xxyz = pbuffer.data(idx_kin_dg + 79);

    auto tk_zz_xxzz = pbuffer.data(idx_kin_dg + 80);

    auto tk_zz_xyyy = pbuffer.data(idx_kin_dg + 81);

    auto tk_zz_xyyz = pbuffer.data(idx_kin_dg + 82);

    auto tk_zz_xyzz = pbuffer.data(idx_kin_dg + 83);

    auto tk_zz_xzzz = pbuffer.data(idx_kin_dg + 84);

    auto tk_zz_yyyy = pbuffer.data(idx_kin_dg + 85);

    auto tk_zz_yyyz = pbuffer.data(idx_kin_dg + 86);

    auto tk_zz_yyzz = pbuffer.data(idx_kin_dg + 87);

    auto tk_zz_yzzz = pbuffer.data(idx_kin_dg + 88);

    auto tk_zz_zzzz = pbuffer.data(idx_kin_dg + 89);

#pragma omp simd aligned(pa_z,           \
                             tk_0_xxxx,  \
                             tk_0_xxxy,  \
                             tk_0_xxxz,  \
                             tk_0_xxyy,  \
                             tk_0_xxyz,  \
                             tk_0_xxzz,  \
                             tk_0_xyyy,  \
                             tk_0_xyyz,  \
                             tk_0_xyzz,  \
                             tk_0_xzzz,  \
                             tk_0_yyyy,  \
                             tk_0_yyyz,  \
                             tk_0_yyzz,  \
                             tk_0_yzzz,  \
                             tk_0_zzzz,  \
                             tk_z_xxx,   \
                             tk_z_xxxx,  \
                             tk_z_xxxy,  \
                             tk_z_xxxz,  \
                             tk_z_xxy,   \
                             tk_z_xxyy,  \
                             tk_z_xxyz,  \
                             tk_z_xxz,   \
                             tk_z_xxzz,  \
                             tk_z_xyy,   \
                             tk_z_xyyy,  \
                             tk_z_xyyz,  \
                             tk_z_xyz,   \
                             tk_z_xyzz,  \
                             tk_z_xzz,   \
                             tk_z_xzzz,  \
                             tk_z_yyy,   \
                             tk_z_yyyy,  \
                             tk_z_yyyz,  \
                             tk_z_yyz,   \
                             tk_z_yyzz,  \
                             tk_z_yzz,   \
                             tk_z_yzzz,  \
                             tk_z_zzz,   \
                             tk_z_zzzz,  \
                             tk_zz_xxxx, \
                             tk_zz_xxxy, \
                             tk_zz_xxxz, \
                             tk_zz_xxyy, \
                             tk_zz_xxyz, \
                             tk_zz_xxzz, \
                             tk_zz_xyyy, \
                             tk_zz_xyyz, \
                             tk_zz_xyzz, \
                             tk_zz_xzzz, \
                             tk_zz_yyyy, \
                             tk_zz_yyyz, \
                             tk_zz_yyzz, \
                             tk_zz_yzzz, \
                             tk_zz_zzzz, \
                             ts_0_xxxx,  \
                             ts_0_xxxy,  \
                             ts_0_xxxz,  \
                             ts_0_xxyy,  \
                             ts_0_xxyz,  \
                             ts_0_xxzz,  \
                             ts_0_xyyy,  \
                             ts_0_xyyz,  \
                             ts_0_xyzz,  \
                             ts_0_xzzz,  \
                             ts_0_yyyy,  \
                             ts_0_yyyz,  \
                             ts_0_yyzz,  \
                             ts_0_yzzz,  \
                             ts_0_zzzz,  \
                             ts_zz_xxxx, \
                             ts_zz_xxxy, \
                             ts_zz_xxxz, \
                             ts_zz_xxyy, \
                             ts_zz_xxyz, \
                             ts_zz_xxzz, \
                             ts_zz_xyyy, \
                             ts_zz_xyyz, \
                             ts_zz_xyzz, \
                             ts_zz_xzzz, \
                             ts_zz_yyyy, \
                             ts_zz_yyyz, \
                             ts_zz_yyzz, \
                             ts_zz_yzzz, \
                             ts_zz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_zz_xxxx[i] = -2.0 * ts_0_xxxx[i] * fbe_0 * fz_0 + tk_0_xxxx[i] * fe_0 + tk_z_xxxx[i] * pa_z[i] + 2.0 * ts_zz_xxxx[i] * fz_0;

        tk_zz_xxxy[i] = -2.0 * ts_0_xxxy[i] * fbe_0 * fz_0 + tk_0_xxxy[i] * fe_0 + tk_z_xxxy[i] * pa_z[i] + 2.0 * ts_zz_xxxy[i] * fz_0;

        tk_zz_xxxz[i] =
            -2.0 * ts_0_xxxz[i] * fbe_0 * fz_0 + tk_0_xxxz[i] * fe_0 + tk_z_xxx[i] * fe_0 + tk_z_xxxz[i] * pa_z[i] + 2.0 * ts_zz_xxxz[i] * fz_0;

        tk_zz_xxyy[i] = -2.0 * ts_0_xxyy[i] * fbe_0 * fz_0 + tk_0_xxyy[i] * fe_0 + tk_z_xxyy[i] * pa_z[i] + 2.0 * ts_zz_xxyy[i] * fz_0;

        tk_zz_xxyz[i] =
            -2.0 * ts_0_xxyz[i] * fbe_0 * fz_0 + tk_0_xxyz[i] * fe_0 + tk_z_xxy[i] * fe_0 + tk_z_xxyz[i] * pa_z[i] + 2.0 * ts_zz_xxyz[i] * fz_0;

        tk_zz_xxzz[i] =
            -2.0 * ts_0_xxzz[i] * fbe_0 * fz_0 + tk_0_xxzz[i] * fe_0 + 2.0 * tk_z_xxz[i] * fe_0 + tk_z_xxzz[i] * pa_z[i] + 2.0 * ts_zz_xxzz[i] * fz_0;

        tk_zz_xyyy[i] = -2.0 * ts_0_xyyy[i] * fbe_0 * fz_0 + tk_0_xyyy[i] * fe_0 + tk_z_xyyy[i] * pa_z[i] + 2.0 * ts_zz_xyyy[i] * fz_0;

        tk_zz_xyyz[i] =
            -2.0 * ts_0_xyyz[i] * fbe_0 * fz_0 + tk_0_xyyz[i] * fe_0 + tk_z_xyy[i] * fe_0 + tk_z_xyyz[i] * pa_z[i] + 2.0 * ts_zz_xyyz[i] * fz_0;

        tk_zz_xyzz[i] =
            -2.0 * ts_0_xyzz[i] * fbe_0 * fz_0 + tk_0_xyzz[i] * fe_0 + 2.0 * tk_z_xyz[i] * fe_0 + tk_z_xyzz[i] * pa_z[i] + 2.0 * ts_zz_xyzz[i] * fz_0;

        tk_zz_xzzz[i] =
            -2.0 * ts_0_xzzz[i] * fbe_0 * fz_0 + tk_0_xzzz[i] * fe_0 + 3.0 * tk_z_xzz[i] * fe_0 + tk_z_xzzz[i] * pa_z[i] + 2.0 * ts_zz_xzzz[i] * fz_0;

        tk_zz_yyyy[i] = -2.0 * ts_0_yyyy[i] * fbe_0 * fz_0 + tk_0_yyyy[i] * fe_0 + tk_z_yyyy[i] * pa_z[i] + 2.0 * ts_zz_yyyy[i] * fz_0;

        tk_zz_yyyz[i] =
            -2.0 * ts_0_yyyz[i] * fbe_0 * fz_0 + tk_0_yyyz[i] * fe_0 + tk_z_yyy[i] * fe_0 + tk_z_yyyz[i] * pa_z[i] + 2.0 * ts_zz_yyyz[i] * fz_0;

        tk_zz_yyzz[i] =
            -2.0 * ts_0_yyzz[i] * fbe_0 * fz_0 + tk_0_yyzz[i] * fe_0 + 2.0 * tk_z_yyz[i] * fe_0 + tk_z_yyzz[i] * pa_z[i] + 2.0 * ts_zz_yyzz[i] * fz_0;

        tk_zz_yzzz[i] =
            -2.0 * ts_0_yzzz[i] * fbe_0 * fz_0 + tk_0_yzzz[i] * fe_0 + 3.0 * tk_z_yzz[i] * fe_0 + tk_z_yzzz[i] * pa_z[i] + 2.0 * ts_zz_yzzz[i] * fz_0;

        tk_zz_zzzz[i] =
            -2.0 * ts_0_zzzz[i] * fbe_0 * fz_0 + tk_0_zzzz[i] * fe_0 + 4.0 * tk_z_zzz[i] * fe_0 + tk_z_zzzz[i] * pa_z[i] + 2.0 * ts_zz_zzzz[i] * fz_0;
    }
}

}  // namespace kinrec
