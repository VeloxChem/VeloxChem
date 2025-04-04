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

#include "ElectricDipoleMomentumPrimRecHP.hpp"

namespace diprec {  // diprec namespace

auto
comp_prim_electric_dipole_momentum_hp(CSimdArray<double>&       pbuffer,
                                      const size_t              idx_dip_hp,
                                      const size_t              idx_dip_fp,
                                      const size_t              idx_dip_gs,
                                      const size_t              idx_ovl_gp,
                                      const size_t              idx_dip_gp,
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

    // Set up components of auxiliary buffer : FP

    auto tr_x_xxx_x = pbuffer.data(idx_dip_fp);

    auto tr_x_xxx_y = pbuffer.data(idx_dip_fp + 1);

    auto tr_x_xxx_z = pbuffer.data(idx_dip_fp + 2);

    auto tr_x_xxy_x = pbuffer.data(idx_dip_fp + 3);

    auto tr_x_xxy_z = pbuffer.data(idx_dip_fp + 5);

    auto tr_x_xxz_x = pbuffer.data(idx_dip_fp + 6);

    auto tr_x_xxz_y = pbuffer.data(idx_dip_fp + 7);

    auto tr_x_xxz_z = pbuffer.data(idx_dip_fp + 8);

    auto tr_x_xyy_x = pbuffer.data(idx_dip_fp + 9);

    auto tr_x_xyy_y = pbuffer.data(idx_dip_fp + 10);

    auto tr_x_xzz_x = pbuffer.data(idx_dip_fp + 15);

    auto tr_x_xzz_z = pbuffer.data(idx_dip_fp + 17);

    auto tr_x_yyy_x = pbuffer.data(idx_dip_fp + 18);

    auto tr_x_yyy_y = pbuffer.data(idx_dip_fp + 19);

    auto tr_x_yyy_z = pbuffer.data(idx_dip_fp + 20);

    auto tr_x_yyz_y = pbuffer.data(idx_dip_fp + 22);

    auto tr_x_yyz_z = pbuffer.data(idx_dip_fp + 23);

    auto tr_x_yzz_x = pbuffer.data(idx_dip_fp + 24);

    auto tr_x_yzz_z = pbuffer.data(idx_dip_fp + 26);

    auto tr_x_zzz_x = pbuffer.data(idx_dip_fp + 27);

    auto tr_x_zzz_y = pbuffer.data(idx_dip_fp + 28);

    auto tr_x_zzz_z = pbuffer.data(idx_dip_fp + 29);

    auto tr_y_xxx_x = pbuffer.data(idx_dip_fp + 30);

    auto tr_y_xxx_y = pbuffer.data(idx_dip_fp + 31);

    auto tr_y_xxx_z = pbuffer.data(idx_dip_fp + 32);

    auto tr_y_xxy_y = pbuffer.data(idx_dip_fp + 34);

    auto tr_y_xxy_z = pbuffer.data(idx_dip_fp + 35);

    auto tr_y_xxz_x = pbuffer.data(idx_dip_fp + 36);

    auto tr_y_xxz_z = pbuffer.data(idx_dip_fp + 38);

    auto tr_y_xyy_x = pbuffer.data(idx_dip_fp + 39);

    auto tr_y_xyy_y = pbuffer.data(idx_dip_fp + 40);

    auto tr_y_xyy_z = pbuffer.data(idx_dip_fp + 41);

    auto tr_y_xyz_z = pbuffer.data(idx_dip_fp + 44);

    auto tr_y_xzz_y = pbuffer.data(idx_dip_fp + 46);

    auto tr_y_xzz_z = pbuffer.data(idx_dip_fp + 47);

    auto tr_y_yyy_x = pbuffer.data(idx_dip_fp + 48);

    auto tr_y_yyy_y = pbuffer.data(idx_dip_fp + 49);

    auto tr_y_yyy_z = pbuffer.data(idx_dip_fp + 50);

    auto tr_y_yyz_x = pbuffer.data(idx_dip_fp + 51);

    auto tr_y_yyz_y = pbuffer.data(idx_dip_fp + 52);

    auto tr_y_yyz_z = pbuffer.data(idx_dip_fp + 53);

    auto tr_y_yzz_y = pbuffer.data(idx_dip_fp + 55);

    auto tr_y_yzz_z = pbuffer.data(idx_dip_fp + 56);

    auto tr_y_zzz_x = pbuffer.data(idx_dip_fp + 57);

    auto tr_y_zzz_y = pbuffer.data(idx_dip_fp + 58);

    auto tr_y_zzz_z = pbuffer.data(idx_dip_fp + 59);

    auto tr_z_xxx_x = pbuffer.data(idx_dip_fp + 60);

    auto tr_z_xxx_y = pbuffer.data(idx_dip_fp + 61);

    auto tr_z_xxx_z = pbuffer.data(idx_dip_fp + 62);

    auto tr_z_xxy_x = pbuffer.data(idx_dip_fp + 63);

    auto tr_z_xxy_y = pbuffer.data(idx_dip_fp + 64);

    auto tr_z_xxz_x = pbuffer.data(idx_dip_fp + 66);

    auto tr_z_xxz_y = pbuffer.data(idx_dip_fp + 67);

    auto tr_z_xxz_z = pbuffer.data(idx_dip_fp + 68);

    auto tr_z_xyy_y = pbuffer.data(idx_dip_fp + 70);

    auto tr_z_xyy_z = pbuffer.data(idx_dip_fp + 71);

    auto tr_z_xyz_y = pbuffer.data(idx_dip_fp + 73);

    auto tr_z_xzz_x = pbuffer.data(idx_dip_fp + 75);

    auto tr_z_xzz_y = pbuffer.data(idx_dip_fp + 76);

    auto tr_z_xzz_z = pbuffer.data(idx_dip_fp + 77);

    auto tr_z_yyy_x = pbuffer.data(idx_dip_fp + 78);

    auto tr_z_yyy_y = pbuffer.data(idx_dip_fp + 79);

    auto tr_z_yyy_z = pbuffer.data(idx_dip_fp + 80);

    auto tr_z_yyz_x = pbuffer.data(idx_dip_fp + 81);

    auto tr_z_yyz_y = pbuffer.data(idx_dip_fp + 82);

    auto tr_z_yyz_z = pbuffer.data(idx_dip_fp + 83);

    auto tr_z_yzz_x = pbuffer.data(idx_dip_fp + 84);

    auto tr_z_yzz_y = pbuffer.data(idx_dip_fp + 85);

    auto tr_z_yzz_z = pbuffer.data(idx_dip_fp + 86);

    auto tr_z_zzz_x = pbuffer.data(idx_dip_fp + 87);

    auto tr_z_zzz_y = pbuffer.data(idx_dip_fp + 88);

    auto tr_z_zzz_z = pbuffer.data(idx_dip_fp + 89);

    // Set up components of auxiliary buffer : GS

    auto tr_x_xxxx_0 = pbuffer.data(idx_dip_gs);

    auto tr_x_xxzz_0 = pbuffer.data(idx_dip_gs + 5);

    auto tr_x_yyyy_0 = pbuffer.data(idx_dip_gs + 10);

    auto tr_x_zzzz_0 = pbuffer.data(idx_dip_gs + 14);

    auto tr_y_xxxx_0 = pbuffer.data(idx_dip_gs + 15);

    auto tr_y_xxyy_0 = pbuffer.data(idx_dip_gs + 18);

    auto tr_y_xyyy_0 = pbuffer.data(idx_dip_gs + 21);

    auto tr_y_yyyy_0 = pbuffer.data(idx_dip_gs + 25);

    auto tr_y_yyzz_0 = pbuffer.data(idx_dip_gs + 27);

    auto tr_y_yzzz_0 = pbuffer.data(idx_dip_gs + 28);

    auto tr_y_zzzz_0 = pbuffer.data(idx_dip_gs + 29);

    auto tr_z_xxxx_0 = pbuffer.data(idx_dip_gs + 30);

    auto tr_z_xxzz_0 = pbuffer.data(idx_dip_gs + 35);

    auto tr_z_xzzz_0 = pbuffer.data(idx_dip_gs + 39);

    auto tr_z_yyyy_0 = pbuffer.data(idx_dip_gs + 40);

    auto tr_z_yyyz_0 = pbuffer.data(idx_dip_gs + 41);

    auto tr_z_yyzz_0 = pbuffer.data(idx_dip_gs + 42);

    auto tr_z_yzzz_0 = pbuffer.data(idx_dip_gs + 43);

    auto tr_z_zzzz_0 = pbuffer.data(idx_dip_gs + 44);

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

    auto ts_yyzz_y = pbuffer.data(idx_ovl_gp + 37);

    auto ts_yyzz_z = pbuffer.data(idx_ovl_gp + 38);

    auto ts_yzzz_y = pbuffer.data(idx_ovl_gp + 40);

    auto ts_yzzz_z = pbuffer.data(idx_ovl_gp + 41);

    auto ts_zzzz_x = pbuffer.data(idx_ovl_gp + 42);

    auto ts_zzzz_y = pbuffer.data(idx_ovl_gp + 43);

    auto ts_zzzz_z = pbuffer.data(idx_ovl_gp + 44);

    // Set up components of auxiliary buffer : GP

    auto tr_x_xxxx_x = pbuffer.data(idx_dip_gp);

    auto tr_x_xxxx_y = pbuffer.data(idx_dip_gp + 1);

    auto tr_x_xxxx_z = pbuffer.data(idx_dip_gp + 2);

    auto tr_x_xxxy_x = pbuffer.data(idx_dip_gp + 3);

    auto tr_x_xxxy_y = pbuffer.data(idx_dip_gp + 4);

    auto tr_x_xxxy_z = pbuffer.data(idx_dip_gp + 5);

    auto tr_x_xxxz_x = pbuffer.data(idx_dip_gp + 6);

    auto tr_x_xxxz_y = pbuffer.data(idx_dip_gp + 7);

    auto tr_x_xxxz_z = pbuffer.data(idx_dip_gp + 8);

    auto tr_x_xxyy_x = pbuffer.data(idx_dip_gp + 9);

    auto tr_x_xxyy_y = pbuffer.data(idx_dip_gp + 10);

    auto tr_x_xxyy_z = pbuffer.data(idx_dip_gp + 11);

    auto tr_x_xxyz_z = pbuffer.data(idx_dip_gp + 14);

    auto tr_x_xxzz_x = pbuffer.data(idx_dip_gp + 15);

    auto tr_x_xxzz_y = pbuffer.data(idx_dip_gp + 16);

    auto tr_x_xxzz_z = pbuffer.data(idx_dip_gp + 17);

    auto tr_x_xyyy_x = pbuffer.data(idx_dip_gp + 18);

    auto tr_x_xyyy_y = pbuffer.data(idx_dip_gp + 19);

    auto tr_x_xyzz_x = pbuffer.data(idx_dip_gp + 24);

    auto tr_x_xzzz_x = pbuffer.data(idx_dip_gp + 27);

    auto tr_x_xzzz_z = pbuffer.data(idx_dip_gp + 29);

    auto tr_x_yyyy_x = pbuffer.data(idx_dip_gp + 30);

    auto tr_x_yyyy_y = pbuffer.data(idx_dip_gp + 31);

    auto tr_x_yyyy_z = pbuffer.data(idx_dip_gp + 32);

    auto tr_x_yyyz_y = pbuffer.data(idx_dip_gp + 34);

    auto tr_x_yyyz_z = pbuffer.data(idx_dip_gp + 35);

    auto tr_x_yyzz_x = pbuffer.data(idx_dip_gp + 36);

    auto tr_x_yyzz_y = pbuffer.data(idx_dip_gp + 37);

    auto tr_x_yyzz_z = pbuffer.data(idx_dip_gp + 38);

    auto tr_x_yzzz_x = pbuffer.data(idx_dip_gp + 39);

    auto tr_x_yzzz_y = pbuffer.data(idx_dip_gp + 40);

    auto tr_x_yzzz_z = pbuffer.data(idx_dip_gp + 41);

    auto tr_x_zzzz_x = pbuffer.data(idx_dip_gp + 42);

    auto tr_x_zzzz_y = pbuffer.data(idx_dip_gp + 43);

    auto tr_x_zzzz_z = pbuffer.data(idx_dip_gp + 44);

    auto tr_y_xxxx_x = pbuffer.data(idx_dip_gp + 45);

    auto tr_y_xxxx_y = pbuffer.data(idx_dip_gp + 46);

    auto tr_y_xxxx_z = pbuffer.data(idx_dip_gp + 47);

    auto tr_y_xxxy_x = pbuffer.data(idx_dip_gp + 48);

    auto tr_y_xxxy_y = pbuffer.data(idx_dip_gp + 49);

    auto tr_y_xxxy_z = pbuffer.data(idx_dip_gp + 50);

    auto tr_y_xxxz_x = pbuffer.data(idx_dip_gp + 51);

    auto tr_y_xxxz_z = pbuffer.data(idx_dip_gp + 53);

    auto tr_y_xxyy_x = pbuffer.data(idx_dip_gp + 54);

    auto tr_y_xxyy_y = pbuffer.data(idx_dip_gp + 55);

    auto tr_y_xxyy_z = pbuffer.data(idx_dip_gp + 56);

    auto tr_y_xxyz_z = pbuffer.data(idx_dip_gp + 59);

    auto tr_y_xxzz_x = pbuffer.data(idx_dip_gp + 60);

    auto tr_y_xxzz_y = pbuffer.data(idx_dip_gp + 61);

    auto tr_y_xxzz_z = pbuffer.data(idx_dip_gp + 62);

    auto tr_y_xyyy_x = pbuffer.data(idx_dip_gp + 63);

    auto tr_y_xyyy_y = pbuffer.data(idx_dip_gp + 64);

    auto tr_y_xyyy_z = pbuffer.data(idx_dip_gp + 65);

    auto tr_y_xyyz_z = pbuffer.data(idx_dip_gp + 68);

    auto tr_y_xyzz_y = pbuffer.data(idx_dip_gp + 70);

    auto tr_y_xyzz_z = pbuffer.data(idx_dip_gp + 71);

    auto tr_y_xzzz_y = pbuffer.data(idx_dip_gp + 73);

    auto tr_y_xzzz_z = pbuffer.data(idx_dip_gp + 74);

    auto tr_y_yyyy_x = pbuffer.data(idx_dip_gp + 75);

    auto tr_y_yyyy_y = pbuffer.data(idx_dip_gp + 76);

    auto tr_y_yyyy_z = pbuffer.data(idx_dip_gp + 77);

    auto tr_y_yyyz_x = pbuffer.data(idx_dip_gp + 78);

    auto tr_y_yyyz_y = pbuffer.data(idx_dip_gp + 79);

    auto tr_y_yyyz_z = pbuffer.data(idx_dip_gp + 80);

    auto tr_y_yyzz_x = pbuffer.data(idx_dip_gp + 81);

    auto tr_y_yyzz_y = pbuffer.data(idx_dip_gp + 82);

    auto tr_y_yyzz_z = pbuffer.data(idx_dip_gp + 83);

    auto tr_y_yzzz_x = pbuffer.data(idx_dip_gp + 84);

    auto tr_y_yzzz_y = pbuffer.data(idx_dip_gp + 85);

    auto tr_y_yzzz_z = pbuffer.data(idx_dip_gp + 86);

    auto tr_y_zzzz_x = pbuffer.data(idx_dip_gp + 87);

    auto tr_y_zzzz_y = pbuffer.data(idx_dip_gp + 88);

    auto tr_y_zzzz_z = pbuffer.data(idx_dip_gp + 89);

    auto tr_z_xxxx_x = pbuffer.data(idx_dip_gp + 90);

    auto tr_z_xxxx_y = pbuffer.data(idx_dip_gp + 91);

    auto tr_z_xxxx_z = pbuffer.data(idx_dip_gp + 92);

    auto tr_z_xxxy_x = pbuffer.data(idx_dip_gp + 93);

    auto tr_z_xxxy_y = pbuffer.data(idx_dip_gp + 94);

    auto tr_z_xxxz_x = pbuffer.data(idx_dip_gp + 96);

    auto tr_z_xxxz_y = pbuffer.data(idx_dip_gp + 97);

    auto tr_z_xxxz_z = pbuffer.data(idx_dip_gp + 98);

    auto tr_z_xxyy_x = pbuffer.data(idx_dip_gp + 99);

    auto tr_z_xxyy_y = pbuffer.data(idx_dip_gp + 100);

    auto tr_z_xxyy_z = pbuffer.data(idx_dip_gp + 101);

    auto tr_z_xxyz_x = pbuffer.data(idx_dip_gp + 102);

    auto tr_z_xxyz_y = pbuffer.data(idx_dip_gp + 103);

    auto tr_z_xxzz_x = pbuffer.data(idx_dip_gp + 105);

    auto tr_z_xxzz_y = pbuffer.data(idx_dip_gp + 106);

    auto tr_z_xxzz_z = pbuffer.data(idx_dip_gp + 107);

    auto tr_z_xyyy_y = pbuffer.data(idx_dip_gp + 109);

    auto tr_z_xyyy_z = pbuffer.data(idx_dip_gp + 110);

    auto tr_z_xyyz_y = pbuffer.data(idx_dip_gp + 112);

    auto tr_z_xyyz_z = pbuffer.data(idx_dip_gp + 113);

    auto tr_z_xyzz_y = pbuffer.data(idx_dip_gp + 115);

    auto tr_z_xzzz_x = pbuffer.data(idx_dip_gp + 117);

    auto tr_z_xzzz_y = pbuffer.data(idx_dip_gp + 118);

    auto tr_z_xzzz_z = pbuffer.data(idx_dip_gp + 119);

    auto tr_z_yyyy_x = pbuffer.data(idx_dip_gp + 120);

    auto tr_z_yyyy_y = pbuffer.data(idx_dip_gp + 121);

    auto tr_z_yyyy_z = pbuffer.data(idx_dip_gp + 122);

    auto tr_z_yyyz_x = pbuffer.data(idx_dip_gp + 123);

    auto tr_z_yyyz_y = pbuffer.data(idx_dip_gp + 124);

    auto tr_z_yyyz_z = pbuffer.data(idx_dip_gp + 125);

    auto tr_z_yyzz_x = pbuffer.data(idx_dip_gp + 126);

    auto tr_z_yyzz_y = pbuffer.data(idx_dip_gp + 127);

    auto tr_z_yyzz_z = pbuffer.data(idx_dip_gp + 128);

    auto tr_z_yzzz_x = pbuffer.data(idx_dip_gp + 129);

    auto tr_z_yzzz_y = pbuffer.data(idx_dip_gp + 130);

    auto tr_z_yzzz_z = pbuffer.data(idx_dip_gp + 131);

    auto tr_z_zzzz_x = pbuffer.data(idx_dip_gp + 132);

    auto tr_z_zzzz_y = pbuffer.data(idx_dip_gp + 133);

    auto tr_z_zzzz_z = pbuffer.data(idx_dip_gp + 134);

    // Set up 0-3 components of targeted buffer : HP

    auto tr_x_xxxxx_x = pbuffer.data(idx_dip_hp);

    auto tr_x_xxxxx_y = pbuffer.data(idx_dip_hp + 1);

    auto tr_x_xxxxx_z = pbuffer.data(idx_dip_hp + 2);

#pragma omp simd aligned(pa_x,             \
                             tr_x_xxx_x,   \
                             tr_x_xxx_y,   \
                             tr_x_xxx_z,   \
                             tr_x_xxxx_0,  \
                             tr_x_xxxx_x,  \
                             tr_x_xxxx_y,  \
                             tr_x_xxxx_z,  \
                             tr_x_xxxxx_x, \
                             tr_x_xxxxx_y, \
                             tr_x_xxxxx_z, \
                             ts_xxxx_x,    \
                             ts_xxxx_y,    \
                             ts_xxxx_z,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxx_x[i] = 4.0 * tr_x_xxx_x[i] * fe_0 + tr_x_xxxx_0[i] * fe_0 + ts_xxxx_x[i] * fe_0 + tr_x_xxxx_x[i] * pa_x[i];

        tr_x_xxxxx_y[i] = 4.0 * tr_x_xxx_y[i] * fe_0 + ts_xxxx_y[i] * fe_0 + tr_x_xxxx_y[i] * pa_x[i];

        tr_x_xxxxx_z[i] = 4.0 * tr_x_xxx_z[i] * fe_0 + ts_xxxx_z[i] * fe_0 + tr_x_xxxx_z[i] * pa_x[i];
    }

    // Set up 3-6 components of targeted buffer : HP

    auto tr_x_xxxxy_x = pbuffer.data(idx_dip_hp + 3);

    auto tr_x_xxxxy_y = pbuffer.data(idx_dip_hp + 4);

    auto tr_x_xxxxy_z = pbuffer.data(idx_dip_hp + 5);

#pragma omp simd aligned(pa_y, tr_x_xxxx_0, tr_x_xxxx_x, tr_x_xxxx_y, tr_x_xxxx_z, tr_x_xxxxy_x, tr_x_xxxxy_y, tr_x_xxxxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxy_x[i] = tr_x_xxxx_x[i] * pa_y[i];

        tr_x_xxxxy_y[i] = tr_x_xxxx_0[i] * fe_0 + tr_x_xxxx_y[i] * pa_y[i];

        tr_x_xxxxy_z[i] = tr_x_xxxx_z[i] * pa_y[i];
    }

    // Set up 6-9 components of targeted buffer : HP

    auto tr_x_xxxxz_x = pbuffer.data(idx_dip_hp + 6);

    auto tr_x_xxxxz_y = pbuffer.data(idx_dip_hp + 7);

    auto tr_x_xxxxz_z = pbuffer.data(idx_dip_hp + 8);

#pragma omp simd aligned(pa_z, tr_x_xxxx_0, tr_x_xxxx_x, tr_x_xxxx_y, tr_x_xxxx_z, tr_x_xxxxz_x, tr_x_xxxxz_y, tr_x_xxxxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxz_x[i] = tr_x_xxxx_x[i] * pa_z[i];

        tr_x_xxxxz_y[i] = tr_x_xxxx_y[i] * pa_z[i];

        tr_x_xxxxz_z[i] = tr_x_xxxx_0[i] * fe_0 + tr_x_xxxx_z[i] * pa_z[i];
    }

    // Set up 9-12 components of targeted buffer : HP

    auto tr_x_xxxyy_x = pbuffer.data(idx_dip_hp + 9);

    auto tr_x_xxxyy_y = pbuffer.data(idx_dip_hp + 10);

    auto tr_x_xxxyy_z = pbuffer.data(idx_dip_hp + 11);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             tr_x_xxx_x,   \
                             tr_x_xxx_z,   \
                             tr_x_xxxy_x,  \
                             tr_x_xxxy_z,  \
                             tr_x_xxxyy_x, \
                             tr_x_xxxyy_y, \
                             tr_x_xxxyy_z, \
                             tr_x_xxyy_y,  \
                             tr_x_xyy_y,   \
                             ts_xxyy_y,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxyy_x[i] = tr_x_xxx_x[i] * fe_0 + tr_x_xxxy_x[i] * pa_y[i];

        tr_x_xxxyy_y[i] = 2.0 * tr_x_xyy_y[i] * fe_0 + ts_xxyy_y[i] * fe_0 + tr_x_xxyy_y[i] * pa_x[i];

        tr_x_xxxyy_z[i] = tr_x_xxx_z[i] * fe_0 + tr_x_xxxy_z[i] * pa_y[i];
    }

    // Set up 12-15 components of targeted buffer : HP

    auto tr_x_xxxyz_x = pbuffer.data(idx_dip_hp + 12);

    auto tr_x_xxxyz_y = pbuffer.data(idx_dip_hp + 13);

    auto tr_x_xxxyz_z = pbuffer.data(idx_dip_hp + 14);

#pragma omp simd aligned(pa_y, pa_z, tr_x_xxxy_y, tr_x_xxxyz_x, tr_x_xxxyz_y, tr_x_xxxyz_z, tr_x_xxxz_x, tr_x_xxxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tr_x_xxxyz_x[i] = tr_x_xxxz_x[i] * pa_y[i];

        tr_x_xxxyz_y[i] = tr_x_xxxy_y[i] * pa_z[i];

        tr_x_xxxyz_z[i] = tr_x_xxxz_z[i] * pa_y[i];
    }

    // Set up 15-18 components of targeted buffer : HP

    auto tr_x_xxxzz_x = pbuffer.data(idx_dip_hp + 15);

    auto tr_x_xxxzz_y = pbuffer.data(idx_dip_hp + 16);

    auto tr_x_xxxzz_z = pbuffer.data(idx_dip_hp + 17);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             tr_x_xxx_x,   \
                             tr_x_xxx_y,   \
                             tr_x_xxxz_x,  \
                             tr_x_xxxz_y,  \
                             tr_x_xxxzz_x, \
                             tr_x_xxxzz_y, \
                             tr_x_xxxzz_z, \
                             tr_x_xxzz_z,  \
                             tr_x_xzz_z,   \
                             ts_xxzz_z,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxzz_x[i] = tr_x_xxx_x[i] * fe_0 + tr_x_xxxz_x[i] * pa_z[i];

        tr_x_xxxzz_y[i] = tr_x_xxx_y[i] * fe_0 + tr_x_xxxz_y[i] * pa_z[i];

        tr_x_xxxzz_z[i] = 2.0 * tr_x_xzz_z[i] * fe_0 + ts_xxzz_z[i] * fe_0 + tr_x_xxzz_z[i] * pa_x[i];
    }

    // Set up 18-21 components of targeted buffer : HP

    auto tr_x_xxyyy_x = pbuffer.data(idx_dip_hp + 18);

    auto tr_x_xxyyy_y = pbuffer.data(idx_dip_hp + 19);

    auto tr_x_xxyyy_z = pbuffer.data(idx_dip_hp + 20);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             tr_x_xxy_x,   \
                             tr_x_xxy_z,   \
                             tr_x_xxyy_x,  \
                             tr_x_xxyy_z,  \
                             tr_x_xxyyy_x, \
                             tr_x_xxyyy_y, \
                             tr_x_xxyyy_z, \
                             tr_x_xyyy_y,  \
                             tr_x_yyy_y,   \
                             ts_xyyy_y,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyyy_x[i] = 2.0 * tr_x_xxy_x[i] * fe_0 + tr_x_xxyy_x[i] * pa_y[i];

        tr_x_xxyyy_y[i] = tr_x_yyy_y[i] * fe_0 + ts_xyyy_y[i] * fe_0 + tr_x_xyyy_y[i] * pa_x[i];

        tr_x_xxyyy_z[i] = 2.0 * tr_x_xxy_z[i] * fe_0 + tr_x_xxyy_z[i] * pa_y[i];
    }

    // Set up 21-24 components of targeted buffer : HP

    auto tr_x_xxyyz_x = pbuffer.data(idx_dip_hp + 21);

    auto tr_x_xxyyz_y = pbuffer.data(idx_dip_hp + 22);

    auto tr_x_xxyyz_z = pbuffer.data(idx_dip_hp + 23);

#pragma omp simd aligned(pa_y, pa_z, tr_x_xxyy_x, tr_x_xxyy_y, tr_x_xxyyz_x, tr_x_xxyyz_y, tr_x_xxyyz_z, tr_x_xxyz_z, tr_x_xxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyyz_x[i] = tr_x_xxyy_x[i] * pa_z[i];

        tr_x_xxyyz_y[i] = tr_x_xxyy_y[i] * pa_z[i];

        tr_x_xxyyz_z[i] = tr_x_xxz_z[i] * fe_0 + tr_x_xxyz_z[i] * pa_y[i];
    }

    // Set up 24-27 components of targeted buffer : HP

    auto tr_x_xxyzz_x = pbuffer.data(idx_dip_hp + 24);

    auto tr_x_xxyzz_y = pbuffer.data(idx_dip_hp + 25);

    auto tr_x_xxyzz_z = pbuffer.data(idx_dip_hp + 26);

#pragma omp simd aligned(pa_y, tr_x_xxyzz_x, tr_x_xxyzz_y, tr_x_xxyzz_z, tr_x_xxzz_0, tr_x_xxzz_x, tr_x_xxzz_y, tr_x_xxzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyzz_x[i] = tr_x_xxzz_x[i] * pa_y[i];

        tr_x_xxyzz_y[i] = tr_x_xxzz_0[i] * fe_0 + tr_x_xxzz_y[i] * pa_y[i];

        tr_x_xxyzz_z[i] = tr_x_xxzz_z[i] * pa_y[i];
    }

    // Set up 27-30 components of targeted buffer : HP

    auto tr_x_xxzzz_x = pbuffer.data(idx_dip_hp + 27);

    auto tr_x_xxzzz_y = pbuffer.data(idx_dip_hp + 28);

    auto tr_x_xxzzz_z = pbuffer.data(idx_dip_hp + 29);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             tr_x_xxz_x,   \
                             tr_x_xxz_y,   \
                             tr_x_xxzz_x,  \
                             tr_x_xxzz_y,  \
                             tr_x_xxzzz_x, \
                             tr_x_xxzzz_y, \
                             tr_x_xxzzz_z, \
                             tr_x_xzzz_z,  \
                             tr_x_zzz_z,   \
                             ts_xzzz_z,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxzzz_x[i] = 2.0 * tr_x_xxz_x[i] * fe_0 + tr_x_xxzz_x[i] * pa_z[i];

        tr_x_xxzzz_y[i] = 2.0 * tr_x_xxz_y[i] * fe_0 + tr_x_xxzz_y[i] * pa_z[i];

        tr_x_xxzzz_z[i] = tr_x_zzz_z[i] * fe_0 + ts_xzzz_z[i] * fe_0 + tr_x_xzzz_z[i] * pa_x[i];
    }

    // Set up 30-33 components of targeted buffer : HP

    auto tr_x_xyyyy_x = pbuffer.data(idx_dip_hp + 30);

    auto tr_x_xyyyy_y = pbuffer.data(idx_dip_hp + 31);

    auto tr_x_xyyyy_z = pbuffer.data(idx_dip_hp + 32);

#pragma omp simd aligned( \
        pa_x, pa_y, tr_x_xyy_x, tr_x_xyyy_x, tr_x_xyyyy_x, tr_x_xyyyy_y, tr_x_xyyyy_z, tr_x_yyyy_y, tr_x_yyyy_z, ts_yyyy_y, ts_yyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyyy_x[i] = 3.0 * tr_x_xyy_x[i] * fe_0 + tr_x_xyyy_x[i] * pa_y[i];

        tr_x_xyyyy_y[i] = ts_yyyy_y[i] * fe_0 + tr_x_yyyy_y[i] * pa_x[i];

        tr_x_xyyyy_z[i] = ts_yyyy_z[i] * fe_0 + tr_x_yyyy_z[i] * pa_x[i];
    }

    // Set up 33-36 components of targeted buffer : HP

    auto tr_x_xyyyz_x = pbuffer.data(idx_dip_hp + 33);

    auto tr_x_xyyyz_y = pbuffer.data(idx_dip_hp + 34);

    auto tr_x_xyyyz_z = pbuffer.data(idx_dip_hp + 35);

#pragma omp simd aligned(pa_x, pa_z, tr_x_xyyy_x, tr_x_xyyy_y, tr_x_xyyyz_x, tr_x_xyyyz_y, tr_x_xyyyz_z, tr_x_yyyz_z, ts_yyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyyz_x[i] = tr_x_xyyy_x[i] * pa_z[i];

        tr_x_xyyyz_y[i] = tr_x_xyyy_y[i] * pa_z[i];

        tr_x_xyyyz_z[i] = ts_yyyz_z[i] * fe_0 + tr_x_yyyz_z[i] * pa_x[i];
    }

    // Set up 36-39 components of targeted buffer : HP

    auto tr_x_xyyzz_x = pbuffer.data(idx_dip_hp + 36);

    auto tr_x_xyyzz_y = pbuffer.data(idx_dip_hp + 37);

    auto tr_x_xyyzz_z = pbuffer.data(idx_dip_hp + 38);

#pragma omp simd aligned( \
        pa_x, pa_y, tr_x_xyyzz_x, tr_x_xyyzz_y, tr_x_xyyzz_z, tr_x_xyzz_x, tr_x_xzz_x, tr_x_yyzz_y, tr_x_yyzz_z, ts_yyzz_y, ts_yyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyzz_x[i] = tr_x_xzz_x[i] * fe_0 + tr_x_xyzz_x[i] * pa_y[i];

        tr_x_xyyzz_y[i] = ts_yyzz_y[i] * fe_0 + tr_x_yyzz_y[i] * pa_x[i];

        tr_x_xyyzz_z[i] = ts_yyzz_z[i] * fe_0 + tr_x_yyzz_z[i] * pa_x[i];
    }

    // Set up 39-42 components of targeted buffer : HP

    auto tr_x_xyzzz_x = pbuffer.data(idx_dip_hp + 39);

    auto tr_x_xyzzz_y = pbuffer.data(idx_dip_hp + 40);

    auto tr_x_xyzzz_z = pbuffer.data(idx_dip_hp + 41);

#pragma omp simd aligned(pa_x, pa_y, tr_x_xyzzz_x, tr_x_xyzzz_y, tr_x_xyzzz_z, tr_x_xzzz_x, tr_x_xzzz_z, tr_x_yzzz_y, ts_yzzz_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyzzz_x[i] = tr_x_xzzz_x[i] * pa_y[i];

        tr_x_xyzzz_y[i] = ts_yzzz_y[i] * fe_0 + tr_x_yzzz_y[i] * pa_x[i];

        tr_x_xyzzz_z[i] = tr_x_xzzz_z[i] * pa_y[i];
    }

    // Set up 42-45 components of targeted buffer : HP

    auto tr_x_xzzzz_x = pbuffer.data(idx_dip_hp + 42);

    auto tr_x_xzzzz_y = pbuffer.data(idx_dip_hp + 43);

    auto tr_x_xzzzz_z = pbuffer.data(idx_dip_hp + 44);

#pragma omp simd aligned( \
        pa_x, pa_z, tr_x_xzz_x, tr_x_xzzz_x, tr_x_xzzzz_x, tr_x_xzzzz_y, tr_x_xzzzz_z, tr_x_zzzz_y, tr_x_zzzz_z, ts_zzzz_y, ts_zzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xzzzz_x[i] = 3.0 * tr_x_xzz_x[i] * fe_0 + tr_x_xzzz_x[i] * pa_z[i];

        tr_x_xzzzz_y[i] = ts_zzzz_y[i] * fe_0 + tr_x_zzzz_y[i] * pa_x[i];

        tr_x_xzzzz_z[i] = ts_zzzz_z[i] * fe_0 + tr_x_zzzz_z[i] * pa_x[i];
    }

    // Set up 45-48 components of targeted buffer : HP

    auto tr_x_yyyyy_x = pbuffer.data(idx_dip_hp + 45);

    auto tr_x_yyyyy_y = pbuffer.data(idx_dip_hp + 46);

    auto tr_x_yyyyy_z = pbuffer.data(idx_dip_hp + 47);

#pragma omp simd aligned(pa_y,             \
                             tr_x_yyy_x,   \
                             tr_x_yyy_y,   \
                             tr_x_yyy_z,   \
                             tr_x_yyyy_0,  \
                             tr_x_yyyy_x,  \
                             tr_x_yyyy_y,  \
                             tr_x_yyyy_z,  \
                             tr_x_yyyyy_x, \
                             tr_x_yyyyy_y, \
                             tr_x_yyyyy_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyyy_x[i] = 4.0 * tr_x_yyy_x[i] * fe_0 + tr_x_yyyy_x[i] * pa_y[i];

        tr_x_yyyyy_y[i] = 4.0 * tr_x_yyy_y[i] * fe_0 + tr_x_yyyy_0[i] * fe_0 + tr_x_yyyy_y[i] * pa_y[i];

        tr_x_yyyyy_z[i] = 4.0 * tr_x_yyy_z[i] * fe_0 + tr_x_yyyy_z[i] * pa_y[i];
    }

    // Set up 48-51 components of targeted buffer : HP

    auto tr_x_yyyyz_x = pbuffer.data(idx_dip_hp + 48);

    auto tr_x_yyyyz_y = pbuffer.data(idx_dip_hp + 49);

    auto tr_x_yyyyz_z = pbuffer.data(idx_dip_hp + 50);

#pragma omp simd aligned(pa_y, pa_z, tr_x_yyyy_x, tr_x_yyyy_y, tr_x_yyyyz_x, tr_x_yyyyz_y, tr_x_yyyyz_z, tr_x_yyyz_z, tr_x_yyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyyz_x[i] = tr_x_yyyy_x[i] * pa_z[i];

        tr_x_yyyyz_y[i] = tr_x_yyyy_y[i] * pa_z[i];

        tr_x_yyyyz_z[i] = 3.0 * tr_x_yyz_z[i] * fe_0 + tr_x_yyyz_z[i] * pa_y[i];
    }

    // Set up 51-54 components of targeted buffer : HP

    auto tr_x_yyyzz_x = pbuffer.data(idx_dip_hp + 51);

    auto tr_x_yyyzz_y = pbuffer.data(idx_dip_hp + 52);

    auto tr_x_yyyzz_z = pbuffer.data(idx_dip_hp + 53);

#pragma omp simd aligned(pa_y,             \
                             pa_z,         \
                             tr_x_yyy_y,   \
                             tr_x_yyyz_y,  \
                             tr_x_yyyzz_x, \
                             tr_x_yyyzz_y, \
                             tr_x_yyyzz_z, \
                             tr_x_yyzz_x,  \
                             tr_x_yyzz_z,  \
                             tr_x_yzz_x,   \
                             tr_x_yzz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyzz_x[i] = 2.0 * tr_x_yzz_x[i] * fe_0 + tr_x_yyzz_x[i] * pa_y[i];

        tr_x_yyyzz_y[i] = tr_x_yyy_y[i] * fe_0 + tr_x_yyyz_y[i] * pa_z[i];

        tr_x_yyyzz_z[i] = 2.0 * tr_x_yzz_z[i] * fe_0 + tr_x_yyzz_z[i] * pa_y[i];
    }

    // Set up 54-57 components of targeted buffer : HP

    auto tr_x_yyzzz_x = pbuffer.data(idx_dip_hp + 54);

    auto tr_x_yyzzz_y = pbuffer.data(idx_dip_hp + 55);

    auto tr_x_yyzzz_z = pbuffer.data(idx_dip_hp + 56);

#pragma omp simd aligned(pa_y,             \
                             pa_z,         \
                             tr_x_yyz_y,   \
                             tr_x_yyzz_y,  \
                             tr_x_yyzzz_x, \
                             tr_x_yyzzz_y, \
                             tr_x_yyzzz_z, \
                             tr_x_yzzz_x,  \
                             tr_x_yzzz_z,  \
                             tr_x_zzz_x,   \
                             tr_x_zzz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyzzz_x[i] = tr_x_zzz_x[i] * fe_0 + tr_x_yzzz_x[i] * pa_y[i];

        tr_x_yyzzz_y[i] = 2.0 * tr_x_yyz_y[i] * fe_0 + tr_x_yyzz_y[i] * pa_z[i];

        tr_x_yyzzz_z[i] = tr_x_zzz_z[i] * fe_0 + tr_x_yzzz_z[i] * pa_y[i];
    }

    // Set up 57-60 components of targeted buffer : HP

    auto tr_x_yzzzz_x = pbuffer.data(idx_dip_hp + 57);

    auto tr_x_yzzzz_y = pbuffer.data(idx_dip_hp + 58);

    auto tr_x_yzzzz_z = pbuffer.data(idx_dip_hp + 59);

#pragma omp simd aligned(pa_y, tr_x_yzzzz_x, tr_x_yzzzz_y, tr_x_yzzzz_z, tr_x_zzzz_0, tr_x_zzzz_x, tr_x_zzzz_y, tr_x_zzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yzzzz_x[i] = tr_x_zzzz_x[i] * pa_y[i];

        tr_x_yzzzz_y[i] = tr_x_zzzz_0[i] * fe_0 + tr_x_zzzz_y[i] * pa_y[i];

        tr_x_yzzzz_z[i] = tr_x_zzzz_z[i] * pa_y[i];
    }

    // Set up 60-63 components of targeted buffer : HP

    auto tr_x_zzzzz_x = pbuffer.data(idx_dip_hp + 60);

    auto tr_x_zzzzz_y = pbuffer.data(idx_dip_hp + 61);

    auto tr_x_zzzzz_z = pbuffer.data(idx_dip_hp + 62);

#pragma omp simd aligned(pa_z,             \
                             tr_x_zzz_x,   \
                             tr_x_zzz_y,   \
                             tr_x_zzz_z,   \
                             tr_x_zzzz_0,  \
                             tr_x_zzzz_x,  \
                             tr_x_zzzz_y,  \
                             tr_x_zzzz_z,  \
                             tr_x_zzzzz_x, \
                             tr_x_zzzzz_y, \
                             tr_x_zzzzz_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_zzzzz_x[i] = 4.0 * tr_x_zzz_x[i] * fe_0 + tr_x_zzzz_x[i] * pa_z[i];

        tr_x_zzzzz_y[i] = 4.0 * tr_x_zzz_y[i] * fe_0 + tr_x_zzzz_y[i] * pa_z[i];

        tr_x_zzzzz_z[i] = 4.0 * tr_x_zzz_z[i] * fe_0 + tr_x_zzzz_0[i] * fe_0 + tr_x_zzzz_z[i] * pa_z[i];
    }

    // Set up 63-66 components of targeted buffer : HP

    auto tr_y_xxxxx_x = pbuffer.data(idx_dip_hp + 63);

    auto tr_y_xxxxx_y = pbuffer.data(idx_dip_hp + 64);

    auto tr_y_xxxxx_z = pbuffer.data(idx_dip_hp + 65);

#pragma omp simd aligned(pa_x,             \
                             tr_y_xxx_x,   \
                             tr_y_xxx_y,   \
                             tr_y_xxx_z,   \
                             tr_y_xxxx_0,  \
                             tr_y_xxxx_x,  \
                             tr_y_xxxx_y,  \
                             tr_y_xxxx_z,  \
                             tr_y_xxxxx_x, \
                             tr_y_xxxxx_y, \
                             tr_y_xxxxx_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxx_x[i] = 4.0 * tr_y_xxx_x[i] * fe_0 + tr_y_xxxx_0[i] * fe_0 + tr_y_xxxx_x[i] * pa_x[i];

        tr_y_xxxxx_y[i] = 4.0 * tr_y_xxx_y[i] * fe_0 + tr_y_xxxx_y[i] * pa_x[i];

        tr_y_xxxxx_z[i] = 4.0 * tr_y_xxx_z[i] * fe_0 + tr_y_xxxx_z[i] * pa_x[i];
    }

    // Set up 66-69 components of targeted buffer : HP

    auto tr_y_xxxxy_x = pbuffer.data(idx_dip_hp + 66);

    auto tr_y_xxxxy_y = pbuffer.data(idx_dip_hp + 67);

    auto tr_y_xxxxy_z = pbuffer.data(idx_dip_hp + 68);

#pragma omp simd aligned( \
        pa_x, pa_y, tr_y_xxxx_x, tr_y_xxxxy_x, tr_y_xxxxy_y, tr_y_xxxxy_z, tr_y_xxxy_y, tr_y_xxxy_z, tr_y_xxy_y, tr_y_xxy_z, ts_xxxx_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxy_x[i] = ts_xxxx_x[i] * fe_0 + tr_y_xxxx_x[i] * pa_y[i];

        tr_y_xxxxy_y[i] = 3.0 * tr_y_xxy_y[i] * fe_0 + tr_y_xxxy_y[i] * pa_x[i];

        tr_y_xxxxy_z[i] = 3.0 * tr_y_xxy_z[i] * fe_0 + tr_y_xxxy_z[i] * pa_x[i];
    }

    // Set up 69-72 components of targeted buffer : HP

    auto tr_y_xxxxz_x = pbuffer.data(idx_dip_hp + 69);

    auto tr_y_xxxxz_y = pbuffer.data(idx_dip_hp + 70);

    auto tr_y_xxxxz_z = pbuffer.data(idx_dip_hp + 71);

#pragma omp simd aligned(pa_x, pa_z, tr_y_xxxx_x, tr_y_xxxx_y, tr_y_xxxxz_x, tr_y_xxxxz_y, tr_y_xxxxz_z, tr_y_xxxz_z, tr_y_xxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxz_x[i] = tr_y_xxxx_x[i] * pa_z[i];

        tr_y_xxxxz_y[i] = tr_y_xxxx_y[i] * pa_z[i];

        tr_y_xxxxz_z[i] = 3.0 * tr_y_xxz_z[i] * fe_0 + tr_y_xxxz_z[i] * pa_x[i];
    }

    // Set up 72-75 components of targeted buffer : HP

    auto tr_y_xxxyy_x = pbuffer.data(idx_dip_hp + 72);

    auto tr_y_xxxyy_y = pbuffer.data(idx_dip_hp + 73);

    auto tr_y_xxxyy_z = pbuffer.data(idx_dip_hp + 74);

#pragma omp simd aligned(pa_x,             \
                             tr_y_xxxyy_x, \
                             tr_y_xxxyy_y, \
                             tr_y_xxxyy_z, \
                             tr_y_xxyy_0,  \
                             tr_y_xxyy_x,  \
                             tr_y_xxyy_y,  \
                             tr_y_xxyy_z,  \
                             tr_y_xyy_x,   \
                             tr_y_xyy_y,   \
                             tr_y_xyy_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxyy_x[i] = 2.0 * tr_y_xyy_x[i] * fe_0 + tr_y_xxyy_0[i] * fe_0 + tr_y_xxyy_x[i] * pa_x[i];

        tr_y_xxxyy_y[i] = 2.0 * tr_y_xyy_y[i] * fe_0 + tr_y_xxyy_y[i] * pa_x[i];

        tr_y_xxxyy_z[i] = 2.0 * tr_y_xyy_z[i] * fe_0 + tr_y_xxyy_z[i] * pa_x[i];
    }

    // Set up 75-78 components of targeted buffer : HP

    auto tr_y_xxxyz_x = pbuffer.data(idx_dip_hp + 75);

    auto tr_y_xxxyz_y = pbuffer.data(idx_dip_hp + 76);

    auto tr_y_xxxyz_z = pbuffer.data(idx_dip_hp + 77);

#pragma omp simd aligned(pa_x, pa_z, tr_y_xxxy_x, tr_y_xxxy_y, tr_y_xxxyz_x, tr_y_xxxyz_y, tr_y_xxxyz_z, tr_y_xxyz_z, tr_y_xyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxyz_x[i] = tr_y_xxxy_x[i] * pa_z[i];

        tr_y_xxxyz_y[i] = tr_y_xxxy_y[i] * pa_z[i];

        tr_y_xxxyz_z[i] = 2.0 * tr_y_xyz_z[i] * fe_0 + tr_y_xxyz_z[i] * pa_x[i];
    }

    // Set up 78-81 components of targeted buffer : HP

    auto tr_y_xxxzz_x = pbuffer.data(idx_dip_hp + 78);

    auto tr_y_xxxzz_y = pbuffer.data(idx_dip_hp + 79);

    auto tr_y_xxxzz_z = pbuffer.data(idx_dip_hp + 80);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             tr_y_xxx_x,   \
                             tr_y_xxxz_x,  \
                             tr_y_xxxzz_x, \
                             tr_y_xxxzz_y, \
                             tr_y_xxxzz_z, \
                             tr_y_xxzz_y,  \
                             tr_y_xxzz_z,  \
                             tr_y_xzz_y,   \
                             tr_y_xzz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxzz_x[i] = tr_y_xxx_x[i] * fe_0 + tr_y_xxxz_x[i] * pa_z[i];

        tr_y_xxxzz_y[i] = 2.0 * tr_y_xzz_y[i] * fe_0 + tr_y_xxzz_y[i] * pa_x[i];

        tr_y_xxxzz_z[i] = 2.0 * tr_y_xzz_z[i] * fe_0 + tr_y_xxzz_z[i] * pa_x[i];
    }

    // Set up 81-84 components of targeted buffer : HP

    auto tr_y_xxyyy_x = pbuffer.data(idx_dip_hp + 81);

    auto tr_y_xxyyy_y = pbuffer.data(idx_dip_hp + 82);

    auto tr_y_xxyyy_z = pbuffer.data(idx_dip_hp + 83);

#pragma omp simd aligned(pa_x,             \
                             tr_y_xxyyy_x, \
                             tr_y_xxyyy_y, \
                             tr_y_xxyyy_z, \
                             tr_y_xyyy_0,  \
                             tr_y_xyyy_x,  \
                             tr_y_xyyy_y,  \
                             tr_y_xyyy_z,  \
                             tr_y_yyy_x,   \
                             tr_y_yyy_y,   \
                             tr_y_yyy_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyyy_x[i] = tr_y_yyy_x[i] * fe_0 + tr_y_xyyy_0[i] * fe_0 + tr_y_xyyy_x[i] * pa_x[i];

        tr_y_xxyyy_y[i] = tr_y_yyy_y[i] * fe_0 + tr_y_xyyy_y[i] * pa_x[i];

        tr_y_xxyyy_z[i] = tr_y_yyy_z[i] * fe_0 + tr_y_xyyy_z[i] * pa_x[i];
    }

    // Set up 84-87 components of targeted buffer : HP

    auto tr_y_xxyyz_x = pbuffer.data(idx_dip_hp + 84);

    auto tr_y_xxyyz_y = pbuffer.data(idx_dip_hp + 85);

    auto tr_y_xxyyz_z = pbuffer.data(idx_dip_hp + 86);

#pragma omp simd aligned(pa_x, pa_z, tr_y_xxyy_x, tr_y_xxyy_y, tr_y_xxyyz_x, tr_y_xxyyz_y, tr_y_xxyyz_z, tr_y_xyyz_z, tr_y_yyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyyz_x[i] = tr_y_xxyy_x[i] * pa_z[i];

        tr_y_xxyyz_y[i] = tr_y_xxyy_y[i] * pa_z[i];

        tr_y_xxyyz_z[i] = tr_y_yyz_z[i] * fe_0 + tr_y_xyyz_z[i] * pa_x[i];
    }

    // Set up 87-90 components of targeted buffer : HP

    auto tr_y_xxyzz_x = pbuffer.data(idx_dip_hp + 87);

    auto tr_y_xxyzz_y = pbuffer.data(idx_dip_hp + 88);

    auto tr_y_xxyzz_z = pbuffer.data(idx_dip_hp + 89);

#pragma omp simd aligned( \
        pa_x, pa_y, tr_y_xxyzz_x, tr_y_xxyzz_y, tr_y_xxyzz_z, tr_y_xxzz_x, tr_y_xyzz_y, tr_y_xyzz_z, tr_y_yzz_y, tr_y_yzz_z, ts_xxzz_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyzz_x[i] = ts_xxzz_x[i] * fe_0 + tr_y_xxzz_x[i] * pa_y[i];

        tr_y_xxyzz_y[i] = tr_y_yzz_y[i] * fe_0 + tr_y_xyzz_y[i] * pa_x[i];

        tr_y_xxyzz_z[i] = tr_y_yzz_z[i] * fe_0 + tr_y_xyzz_z[i] * pa_x[i];
    }

    // Set up 90-93 components of targeted buffer : HP

    auto tr_y_xxzzz_x = pbuffer.data(idx_dip_hp + 90);

    auto tr_y_xxzzz_y = pbuffer.data(idx_dip_hp + 91);

    auto tr_y_xxzzz_z = pbuffer.data(idx_dip_hp + 92);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             tr_y_xxz_x,   \
                             tr_y_xxzz_x,  \
                             tr_y_xxzzz_x, \
                             tr_y_xxzzz_y, \
                             tr_y_xxzzz_z, \
                             tr_y_xzzz_y,  \
                             tr_y_xzzz_z,  \
                             tr_y_zzz_y,   \
                             tr_y_zzz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxzzz_x[i] = 2.0 * tr_y_xxz_x[i] * fe_0 + tr_y_xxzz_x[i] * pa_z[i];

        tr_y_xxzzz_y[i] = tr_y_zzz_y[i] * fe_0 + tr_y_xzzz_y[i] * pa_x[i];

        tr_y_xxzzz_z[i] = tr_y_zzz_z[i] * fe_0 + tr_y_xzzz_z[i] * pa_x[i];
    }

    // Set up 93-96 components of targeted buffer : HP

    auto tr_y_xyyyy_x = pbuffer.data(idx_dip_hp + 93);

    auto tr_y_xyyyy_y = pbuffer.data(idx_dip_hp + 94);

    auto tr_y_xyyyy_z = pbuffer.data(idx_dip_hp + 95);

#pragma omp simd aligned(pa_x, tr_y_xyyyy_x, tr_y_xyyyy_y, tr_y_xyyyy_z, tr_y_yyyy_0, tr_y_yyyy_x, tr_y_yyyy_y, tr_y_yyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyyy_x[i] = tr_y_yyyy_0[i] * fe_0 + tr_y_yyyy_x[i] * pa_x[i];

        tr_y_xyyyy_y[i] = tr_y_yyyy_y[i] * pa_x[i];

        tr_y_xyyyy_z[i] = tr_y_yyyy_z[i] * pa_x[i];
    }

    // Set up 96-99 components of targeted buffer : HP

    auto tr_y_xyyyz_x = pbuffer.data(idx_dip_hp + 96);

    auto tr_y_xyyyz_y = pbuffer.data(idx_dip_hp + 97);

    auto tr_y_xyyyz_z = pbuffer.data(idx_dip_hp + 98);

#pragma omp simd aligned(pa_x, pa_z, tr_y_xyyy_x, tr_y_xyyyz_x, tr_y_xyyyz_y, tr_y_xyyyz_z, tr_y_yyyz_y, tr_y_yyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tr_y_xyyyz_x[i] = tr_y_xyyy_x[i] * pa_z[i];

        tr_y_xyyyz_y[i] = tr_y_yyyz_y[i] * pa_x[i];

        tr_y_xyyyz_z[i] = tr_y_yyyz_z[i] * pa_x[i];
    }

    // Set up 99-102 components of targeted buffer : HP

    auto tr_y_xyyzz_x = pbuffer.data(idx_dip_hp + 99);

    auto tr_y_xyyzz_y = pbuffer.data(idx_dip_hp + 100);

    auto tr_y_xyyzz_z = pbuffer.data(idx_dip_hp + 101);

#pragma omp simd aligned(pa_x, tr_y_xyyzz_x, tr_y_xyyzz_y, tr_y_xyyzz_z, tr_y_yyzz_0, tr_y_yyzz_x, tr_y_yyzz_y, tr_y_yyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyzz_x[i] = tr_y_yyzz_0[i] * fe_0 + tr_y_yyzz_x[i] * pa_x[i];

        tr_y_xyyzz_y[i] = tr_y_yyzz_y[i] * pa_x[i];

        tr_y_xyyzz_z[i] = tr_y_yyzz_z[i] * pa_x[i];
    }

    // Set up 102-105 components of targeted buffer : HP

    auto tr_y_xyzzz_x = pbuffer.data(idx_dip_hp + 102);

    auto tr_y_xyzzz_y = pbuffer.data(idx_dip_hp + 103);

    auto tr_y_xyzzz_z = pbuffer.data(idx_dip_hp + 104);

#pragma omp simd aligned(pa_x, tr_y_xyzzz_x, tr_y_xyzzz_y, tr_y_xyzzz_z, tr_y_yzzz_0, tr_y_yzzz_x, tr_y_yzzz_y, tr_y_yzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyzzz_x[i] = tr_y_yzzz_0[i] * fe_0 + tr_y_yzzz_x[i] * pa_x[i];

        tr_y_xyzzz_y[i] = tr_y_yzzz_y[i] * pa_x[i];

        tr_y_xyzzz_z[i] = tr_y_yzzz_z[i] * pa_x[i];
    }

    // Set up 105-108 components of targeted buffer : HP

    auto tr_y_xzzzz_x = pbuffer.data(idx_dip_hp + 105);

    auto tr_y_xzzzz_y = pbuffer.data(idx_dip_hp + 106);

    auto tr_y_xzzzz_z = pbuffer.data(idx_dip_hp + 107);

#pragma omp simd aligned(pa_x, tr_y_xzzzz_x, tr_y_xzzzz_y, tr_y_xzzzz_z, tr_y_zzzz_0, tr_y_zzzz_x, tr_y_zzzz_y, tr_y_zzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xzzzz_x[i] = tr_y_zzzz_0[i] * fe_0 + tr_y_zzzz_x[i] * pa_x[i];

        tr_y_xzzzz_y[i] = tr_y_zzzz_y[i] * pa_x[i];

        tr_y_xzzzz_z[i] = tr_y_zzzz_z[i] * pa_x[i];
    }

    // Set up 108-111 components of targeted buffer : HP

    auto tr_y_yyyyy_x = pbuffer.data(idx_dip_hp + 108);

    auto tr_y_yyyyy_y = pbuffer.data(idx_dip_hp + 109);

    auto tr_y_yyyyy_z = pbuffer.data(idx_dip_hp + 110);

#pragma omp simd aligned(pa_y,             \
                             tr_y_yyy_x,   \
                             tr_y_yyy_y,   \
                             tr_y_yyy_z,   \
                             tr_y_yyyy_0,  \
                             tr_y_yyyy_x,  \
                             tr_y_yyyy_y,  \
                             tr_y_yyyy_z,  \
                             tr_y_yyyyy_x, \
                             tr_y_yyyyy_y, \
                             tr_y_yyyyy_z, \
                             ts_yyyy_x,    \
                             ts_yyyy_y,    \
                             ts_yyyy_z,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyyy_x[i] = 4.0 * tr_y_yyy_x[i] * fe_0 + ts_yyyy_x[i] * fe_0 + tr_y_yyyy_x[i] * pa_y[i];

        tr_y_yyyyy_y[i] = 4.0 * tr_y_yyy_y[i] * fe_0 + tr_y_yyyy_0[i] * fe_0 + ts_yyyy_y[i] * fe_0 + tr_y_yyyy_y[i] * pa_y[i];

        tr_y_yyyyy_z[i] = 4.0 * tr_y_yyy_z[i] * fe_0 + ts_yyyy_z[i] * fe_0 + tr_y_yyyy_z[i] * pa_y[i];
    }

    // Set up 111-114 components of targeted buffer : HP

    auto tr_y_yyyyz_x = pbuffer.data(idx_dip_hp + 111);

    auto tr_y_yyyyz_y = pbuffer.data(idx_dip_hp + 112);

    auto tr_y_yyyyz_z = pbuffer.data(idx_dip_hp + 113);

#pragma omp simd aligned(pa_z, tr_y_yyyy_0, tr_y_yyyy_x, tr_y_yyyy_y, tr_y_yyyy_z, tr_y_yyyyz_x, tr_y_yyyyz_y, tr_y_yyyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyyz_x[i] = tr_y_yyyy_x[i] * pa_z[i];

        tr_y_yyyyz_y[i] = tr_y_yyyy_y[i] * pa_z[i];

        tr_y_yyyyz_z[i] = tr_y_yyyy_0[i] * fe_0 + tr_y_yyyy_z[i] * pa_z[i];
    }

    // Set up 114-117 components of targeted buffer : HP

    auto tr_y_yyyzz_x = pbuffer.data(idx_dip_hp + 114);

    auto tr_y_yyyzz_y = pbuffer.data(idx_dip_hp + 115);

    auto tr_y_yyyzz_z = pbuffer.data(idx_dip_hp + 116);

#pragma omp simd aligned(pa_y,             \
                             pa_z,         \
                             tr_y_yyy_x,   \
                             tr_y_yyy_y,   \
                             tr_y_yyyz_x,  \
                             tr_y_yyyz_y,  \
                             tr_y_yyyzz_x, \
                             tr_y_yyyzz_y, \
                             tr_y_yyyzz_z, \
                             tr_y_yyzz_z,  \
                             tr_y_yzz_z,   \
                             ts_yyzz_z,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyzz_x[i] = tr_y_yyy_x[i] * fe_0 + tr_y_yyyz_x[i] * pa_z[i];

        tr_y_yyyzz_y[i] = tr_y_yyy_y[i] * fe_0 + tr_y_yyyz_y[i] * pa_z[i];

        tr_y_yyyzz_z[i] = 2.0 * tr_y_yzz_z[i] * fe_0 + ts_yyzz_z[i] * fe_0 + tr_y_yyzz_z[i] * pa_y[i];
    }

    // Set up 117-120 components of targeted buffer : HP

    auto tr_y_yyzzz_x = pbuffer.data(idx_dip_hp + 117);

    auto tr_y_yyzzz_y = pbuffer.data(idx_dip_hp + 118);

    auto tr_y_yyzzz_z = pbuffer.data(idx_dip_hp + 119);

#pragma omp simd aligned(pa_y,             \
                             pa_z,         \
                             tr_y_yyz_x,   \
                             tr_y_yyz_y,   \
                             tr_y_yyzz_x,  \
                             tr_y_yyzz_y,  \
                             tr_y_yyzzz_x, \
                             tr_y_yyzzz_y, \
                             tr_y_yyzzz_z, \
                             tr_y_yzzz_z,  \
                             tr_y_zzz_z,   \
                             ts_yzzz_z,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyzzz_x[i] = 2.0 * tr_y_yyz_x[i] * fe_0 + tr_y_yyzz_x[i] * pa_z[i];

        tr_y_yyzzz_y[i] = 2.0 * tr_y_yyz_y[i] * fe_0 + tr_y_yyzz_y[i] * pa_z[i];

        tr_y_yyzzz_z[i] = tr_y_zzz_z[i] * fe_0 + ts_yzzz_z[i] * fe_0 + tr_y_yzzz_z[i] * pa_y[i];
    }

    // Set up 120-123 components of targeted buffer : HP

    auto tr_y_yzzzz_x = pbuffer.data(idx_dip_hp + 120);

    auto tr_y_yzzzz_y = pbuffer.data(idx_dip_hp + 121);

    auto tr_y_yzzzz_z = pbuffer.data(idx_dip_hp + 122);

#pragma omp simd aligned( \
        pa_y, pa_z, tr_y_yzz_y, tr_y_yzzz_y, tr_y_yzzzz_x, tr_y_yzzzz_y, tr_y_yzzzz_z, tr_y_zzzz_x, tr_y_zzzz_z, ts_zzzz_x, ts_zzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yzzzz_x[i] = ts_zzzz_x[i] * fe_0 + tr_y_zzzz_x[i] * pa_y[i];

        tr_y_yzzzz_y[i] = 3.0 * tr_y_yzz_y[i] * fe_0 + tr_y_yzzz_y[i] * pa_z[i];

        tr_y_yzzzz_z[i] = ts_zzzz_z[i] * fe_0 + tr_y_zzzz_z[i] * pa_y[i];
    }

    // Set up 123-126 components of targeted buffer : HP

    auto tr_y_zzzzz_x = pbuffer.data(idx_dip_hp + 123);

    auto tr_y_zzzzz_y = pbuffer.data(idx_dip_hp + 124);

    auto tr_y_zzzzz_z = pbuffer.data(idx_dip_hp + 125);

#pragma omp simd aligned(pa_z,             \
                             tr_y_zzz_x,   \
                             tr_y_zzz_y,   \
                             tr_y_zzz_z,   \
                             tr_y_zzzz_0,  \
                             tr_y_zzzz_x,  \
                             tr_y_zzzz_y,  \
                             tr_y_zzzz_z,  \
                             tr_y_zzzzz_x, \
                             tr_y_zzzzz_y, \
                             tr_y_zzzzz_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_zzzzz_x[i] = 4.0 * tr_y_zzz_x[i] * fe_0 + tr_y_zzzz_x[i] * pa_z[i];

        tr_y_zzzzz_y[i] = 4.0 * tr_y_zzz_y[i] * fe_0 + tr_y_zzzz_y[i] * pa_z[i];

        tr_y_zzzzz_z[i] = 4.0 * tr_y_zzz_z[i] * fe_0 + tr_y_zzzz_0[i] * fe_0 + tr_y_zzzz_z[i] * pa_z[i];
    }

    // Set up 126-129 components of targeted buffer : HP

    auto tr_z_xxxxx_x = pbuffer.data(idx_dip_hp + 126);

    auto tr_z_xxxxx_y = pbuffer.data(idx_dip_hp + 127);

    auto tr_z_xxxxx_z = pbuffer.data(idx_dip_hp + 128);

#pragma omp simd aligned(pa_x,             \
                             tr_z_xxx_x,   \
                             tr_z_xxx_y,   \
                             tr_z_xxx_z,   \
                             tr_z_xxxx_0,  \
                             tr_z_xxxx_x,  \
                             tr_z_xxxx_y,  \
                             tr_z_xxxx_z,  \
                             tr_z_xxxxx_x, \
                             tr_z_xxxxx_y, \
                             tr_z_xxxxx_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxx_x[i] = 4.0 * tr_z_xxx_x[i] * fe_0 + tr_z_xxxx_0[i] * fe_0 + tr_z_xxxx_x[i] * pa_x[i];

        tr_z_xxxxx_y[i] = 4.0 * tr_z_xxx_y[i] * fe_0 + tr_z_xxxx_y[i] * pa_x[i];

        tr_z_xxxxx_z[i] = 4.0 * tr_z_xxx_z[i] * fe_0 + tr_z_xxxx_z[i] * pa_x[i];
    }

    // Set up 129-132 components of targeted buffer : HP

    auto tr_z_xxxxy_x = pbuffer.data(idx_dip_hp + 129);

    auto tr_z_xxxxy_y = pbuffer.data(idx_dip_hp + 130);

    auto tr_z_xxxxy_z = pbuffer.data(idx_dip_hp + 131);

#pragma omp simd aligned(pa_x, pa_y, tr_z_xxxx_x, tr_z_xxxx_z, tr_z_xxxxy_x, tr_z_xxxxy_y, tr_z_xxxxy_z, tr_z_xxxy_y, tr_z_xxy_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxy_x[i] = tr_z_xxxx_x[i] * pa_y[i];

        tr_z_xxxxy_y[i] = 3.0 * tr_z_xxy_y[i] * fe_0 + tr_z_xxxy_y[i] * pa_x[i];

        tr_z_xxxxy_z[i] = tr_z_xxxx_z[i] * pa_y[i];
    }

    // Set up 132-135 components of targeted buffer : HP

    auto tr_z_xxxxz_x = pbuffer.data(idx_dip_hp + 132);

    auto tr_z_xxxxz_y = pbuffer.data(idx_dip_hp + 133);

    auto tr_z_xxxxz_z = pbuffer.data(idx_dip_hp + 134);

#pragma omp simd aligned( \
        pa_x, pa_z, tr_z_xxxx_x, tr_z_xxxxz_x, tr_z_xxxxz_y, tr_z_xxxxz_z, tr_z_xxxz_y, tr_z_xxxz_z, tr_z_xxz_y, tr_z_xxz_z, ts_xxxx_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxz_x[i] = ts_xxxx_x[i] * fe_0 + tr_z_xxxx_x[i] * pa_z[i];

        tr_z_xxxxz_y[i] = 3.0 * tr_z_xxz_y[i] * fe_0 + tr_z_xxxz_y[i] * pa_x[i];

        tr_z_xxxxz_z[i] = 3.0 * tr_z_xxz_z[i] * fe_0 + tr_z_xxxz_z[i] * pa_x[i];
    }

    // Set up 135-138 components of targeted buffer : HP

    auto tr_z_xxxyy_x = pbuffer.data(idx_dip_hp + 135);

    auto tr_z_xxxyy_y = pbuffer.data(idx_dip_hp + 136);

    auto tr_z_xxxyy_z = pbuffer.data(idx_dip_hp + 137);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             tr_z_xxx_x,   \
                             tr_z_xxxy_x,  \
                             tr_z_xxxyy_x, \
                             tr_z_xxxyy_y, \
                             tr_z_xxxyy_z, \
                             tr_z_xxyy_y,  \
                             tr_z_xxyy_z,  \
                             tr_z_xyy_y,   \
                             tr_z_xyy_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxyy_x[i] = tr_z_xxx_x[i] * fe_0 + tr_z_xxxy_x[i] * pa_y[i];

        tr_z_xxxyy_y[i] = 2.0 * tr_z_xyy_y[i] * fe_0 + tr_z_xxyy_y[i] * pa_x[i];

        tr_z_xxxyy_z[i] = 2.0 * tr_z_xyy_z[i] * fe_0 + tr_z_xxyy_z[i] * pa_x[i];
    }

    // Set up 138-141 components of targeted buffer : HP

    auto tr_z_xxxyz_x = pbuffer.data(idx_dip_hp + 138);

    auto tr_z_xxxyz_y = pbuffer.data(idx_dip_hp + 139);

    auto tr_z_xxxyz_z = pbuffer.data(idx_dip_hp + 140);

#pragma omp simd aligned(pa_x, pa_y, tr_z_xxxyz_x, tr_z_xxxyz_y, tr_z_xxxyz_z, tr_z_xxxz_x, tr_z_xxxz_z, tr_z_xxyz_y, tr_z_xyz_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxyz_x[i] = tr_z_xxxz_x[i] * pa_y[i];

        tr_z_xxxyz_y[i] = 2.0 * tr_z_xyz_y[i] * fe_0 + tr_z_xxyz_y[i] * pa_x[i];

        tr_z_xxxyz_z[i] = tr_z_xxxz_z[i] * pa_y[i];
    }

    // Set up 141-144 components of targeted buffer : HP

    auto tr_z_xxxzz_x = pbuffer.data(idx_dip_hp + 141);

    auto tr_z_xxxzz_y = pbuffer.data(idx_dip_hp + 142);

    auto tr_z_xxxzz_z = pbuffer.data(idx_dip_hp + 143);

#pragma omp simd aligned(pa_x,             \
                             tr_z_xxxzz_x, \
                             tr_z_xxxzz_y, \
                             tr_z_xxxzz_z, \
                             tr_z_xxzz_0,  \
                             tr_z_xxzz_x,  \
                             tr_z_xxzz_y,  \
                             tr_z_xxzz_z,  \
                             tr_z_xzz_x,   \
                             tr_z_xzz_y,   \
                             tr_z_xzz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxzz_x[i] = 2.0 * tr_z_xzz_x[i] * fe_0 + tr_z_xxzz_0[i] * fe_0 + tr_z_xxzz_x[i] * pa_x[i];

        tr_z_xxxzz_y[i] = 2.0 * tr_z_xzz_y[i] * fe_0 + tr_z_xxzz_y[i] * pa_x[i];

        tr_z_xxxzz_z[i] = 2.0 * tr_z_xzz_z[i] * fe_0 + tr_z_xxzz_z[i] * pa_x[i];
    }

    // Set up 144-147 components of targeted buffer : HP

    auto tr_z_xxyyy_x = pbuffer.data(idx_dip_hp + 144);

    auto tr_z_xxyyy_y = pbuffer.data(idx_dip_hp + 145);

    auto tr_z_xxyyy_z = pbuffer.data(idx_dip_hp + 146);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             tr_z_xxy_x,   \
                             tr_z_xxyy_x,  \
                             tr_z_xxyyy_x, \
                             tr_z_xxyyy_y, \
                             tr_z_xxyyy_z, \
                             tr_z_xyyy_y,  \
                             tr_z_xyyy_z,  \
                             tr_z_yyy_y,   \
                             tr_z_yyy_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyyy_x[i] = 2.0 * tr_z_xxy_x[i] * fe_0 + tr_z_xxyy_x[i] * pa_y[i];

        tr_z_xxyyy_y[i] = tr_z_yyy_y[i] * fe_0 + tr_z_xyyy_y[i] * pa_x[i];

        tr_z_xxyyy_z[i] = tr_z_yyy_z[i] * fe_0 + tr_z_xyyy_z[i] * pa_x[i];
    }

    // Set up 147-150 components of targeted buffer : HP

    auto tr_z_xxyyz_x = pbuffer.data(idx_dip_hp + 147);

    auto tr_z_xxyyz_y = pbuffer.data(idx_dip_hp + 148);

    auto tr_z_xxyyz_z = pbuffer.data(idx_dip_hp + 149);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             tr_z_xxyyz_x, \
                             tr_z_xxyyz_y, \
                             tr_z_xxyyz_z, \
                             tr_z_xxyz_x,  \
                             tr_z_xxz_x,   \
                             tr_z_xyyz_y,  \
                             tr_z_xyyz_z,  \
                             tr_z_yyz_y,   \
                             tr_z_yyz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyyz_x[i] = tr_z_xxz_x[i] * fe_0 + tr_z_xxyz_x[i] * pa_y[i];

        tr_z_xxyyz_y[i] = tr_z_yyz_y[i] * fe_0 + tr_z_xyyz_y[i] * pa_x[i];

        tr_z_xxyyz_z[i] = tr_z_yyz_z[i] * fe_0 + tr_z_xyyz_z[i] * pa_x[i];
    }

    // Set up 150-153 components of targeted buffer : HP

    auto tr_z_xxyzz_x = pbuffer.data(idx_dip_hp + 150);

    auto tr_z_xxyzz_y = pbuffer.data(idx_dip_hp + 151);

    auto tr_z_xxyzz_z = pbuffer.data(idx_dip_hp + 152);

#pragma omp simd aligned(pa_x, pa_y, tr_z_xxyzz_x, tr_z_xxyzz_y, tr_z_xxyzz_z, tr_z_xxzz_x, tr_z_xxzz_z, tr_z_xyzz_y, tr_z_yzz_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyzz_x[i] = tr_z_xxzz_x[i] * pa_y[i];

        tr_z_xxyzz_y[i] = tr_z_yzz_y[i] * fe_0 + tr_z_xyzz_y[i] * pa_x[i];

        tr_z_xxyzz_z[i] = tr_z_xxzz_z[i] * pa_y[i];
    }

    // Set up 153-156 components of targeted buffer : HP

    auto tr_z_xxzzz_x = pbuffer.data(idx_dip_hp + 153);

    auto tr_z_xxzzz_y = pbuffer.data(idx_dip_hp + 154);

    auto tr_z_xxzzz_z = pbuffer.data(idx_dip_hp + 155);

#pragma omp simd aligned(pa_x,             \
                             tr_z_xxzzz_x, \
                             tr_z_xxzzz_y, \
                             tr_z_xxzzz_z, \
                             tr_z_xzzz_0,  \
                             tr_z_xzzz_x,  \
                             tr_z_xzzz_y,  \
                             tr_z_xzzz_z,  \
                             tr_z_zzz_x,   \
                             tr_z_zzz_y,   \
                             tr_z_zzz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxzzz_x[i] = tr_z_zzz_x[i] * fe_0 + tr_z_xzzz_0[i] * fe_0 + tr_z_xzzz_x[i] * pa_x[i];

        tr_z_xxzzz_y[i] = tr_z_zzz_y[i] * fe_0 + tr_z_xzzz_y[i] * pa_x[i];

        tr_z_xxzzz_z[i] = tr_z_zzz_z[i] * fe_0 + tr_z_xzzz_z[i] * pa_x[i];
    }

    // Set up 156-159 components of targeted buffer : HP

    auto tr_z_xyyyy_x = pbuffer.data(idx_dip_hp + 156);

    auto tr_z_xyyyy_y = pbuffer.data(idx_dip_hp + 157);

    auto tr_z_xyyyy_z = pbuffer.data(idx_dip_hp + 158);

#pragma omp simd aligned(pa_x, tr_z_xyyyy_x, tr_z_xyyyy_y, tr_z_xyyyy_z, tr_z_yyyy_0, tr_z_yyyy_x, tr_z_yyyy_y, tr_z_yyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyyy_x[i] = tr_z_yyyy_0[i] * fe_0 + tr_z_yyyy_x[i] * pa_x[i];

        tr_z_xyyyy_y[i] = tr_z_yyyy_y[i] * pa_x[i];

        tr_z_xyyyy_z[i] = tr_z_yyyy_z[i] * pa_x[i];
    }

    // Set up 159-162 components of targeted buffer : HP

    auto tr_z_xyyyz_x = pbuffer.data(idx_dip_hp + 159);

    auto tr_z_xyyyz_y = pbuffer.data(idx_dip_hp + 160);

    auto tr_z_xyyyz_z = pbuffer.data(idx_dip_hp + 161);

#pragma omp simd aligned(pa_x, tr_z_xyyyz_x, tr_z_xyyyz_y, tr_z_xyyyz_z, tr_z_yyyz_0, tr_z_yyyz_x, tr_z_yyyz_y, tr_z_yyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyyz_x[i] = tr_z_yyyz_0[i] * fe_0 + tr_z_yyyz_x[i] * pa_x[i];

        tr_z_xyyyz_y[i] = tr_z_yyyz_y[i] * pa_x[i];

        tr_z_xyyyz_z[i] = tr_z_yyyz_z[i] * pa_x[i];
    }

    // Set up 162-165 components of targeted buffer : HP

    auto tr_z_xyyzz_x = pbuffer.data(idx_dip_hp + 162);

    auto tr_z_xyyzz_y = pbuffer.data(idx_dip_hp + 163);

    auto tr_z_xyyzz_z = pbuffer.data(idx_dip_hp + 164);

#pragma omp simd aligned(pa_x, tr_z_xyyzz_x, tr_z_xyyzz_y, tr_z_xyyzz_z, tr_z_yyzz_0, tr_z_yyzz_x, tr_z_yyzz_y, tr_z_yyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyzz_x[i] = tr_z_yyzz_0[i] * fe_0 + tr_z_yyzz_x[i] * pa_x[i];

        tr_z_xyyzz_y[i] = tr_z_yyzz_y[i] * pa_x[i];

        tr_z_xyyzz_z[i] = tr_z_yyzz_z[i] * pa_x[i];
    }

    // Set up 165-168 components of targeted buffer : HP

    auto tr_z_xyzzz_x = pbuffer.data(idx_dip_hp + 165);

    auto tr_z_xyzzz_y = pbuffer.data(idx_dip_hp + 166);

    auto tr_z_xyzzz_z = pbuffer.data(idx_dip_hp + 167);

#pragma omp simd aligned(pa_x, pa_y, tr_z_xyzzz_x, tr_z_xyzzz_y, tr_z_xyzzz_z, tr_z_xzzz_x, tr_z_yzzz_y, tr_z_yzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tr_z_xyzzz_x[i] = tr_z_xzzz_x[i] * pa_y[i];

        tr_z_xyzzz_y[i] = tr_z_yzzz_y[i] * pa_x[i];

        tr_z_xyzzz_z[i] = tr_z_yzzz_z[i] * pa_x[i];
    }

    // Set up 168-171 components of targeted buffer : HP

    auto tr_z_xzzzz_x = pbuffer.data(idx_dip_hp + 168);

    auto tr_z_xzzzz_y = pbuffer.data(idx_dip_hp + 169);

    auto tr_z_xzzzz_z = pbuffer.data(idx_dip_hp + 170);

#pragma omp simd aligned(pa_x, tr_z_xzzzz_x, tr_z_xzzzz_y, tr_z_xzzzz_z, tr_z_zzzz_0, tr_z_zzzz_x, tr_z_zzzz_y, tr_z_zzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xzzzz_x[i] = tr_z_zzzz_0[i] * fe_0 + tr_z_zzzz_x[i] * pa_x[i];

        tr_z_xzzzz_y[i] = tr_z_zzzz_y[i] * pa_x[i];

        tr_z_xzzzz_z[i] = tr_z_zzzz_z[i] * pa_x[i];
    }

    // Set up 171-174 components of targeted buffer : HP

    auto tr_z_yyyyy_x = pbuffer.data(idx_dip_hp + 171);

    auto tr_z_yyyyy_y = pbuffer.data(idx_dip_hp + 172);

    auto tr_z_yyyyy_z = pbuffer.data(idx_dip_hp + 173);

#pragma omp simd aligned(pa_y,             \
                             tr_z_yyy_x,   \
                             tr_z_yyy_y,   \
                             tr_z_yyy_z,   \
                             tr_z_yyyy_0,  \
                             tr_z_yyyy_x,  \
                             tr_z_yyyy_y,  \
                             tr_z_yyyy_z,  \
                             tr_z_yyyyy_x, \
                             tr_z_yyyyy_y, \
                             tr_z_yyyyy_z, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyyy_x[i] = 4.0 * tr_z_yyy_x[i] * fe_0 + tr_z_yyyy_x[i] * pa_y[i];

        tr_z_yyyyy_y[i] = 4.0 * tr_z_yyy_y[i] * fe_0 + tr_z_yyyy_0[i] * fe_0 + tr_z_yyyy_y[i] * pa_y[i];

        tr_z_yyyyy_z[i] = 4.0 * tr_z_yyy_z[i] * fe_0 + tr_z_yyyy_z[i] * pa_y[i];
    }

    // Set up 174-177 components of targeted buffer : HP

    auto tr_z_yyyyz_x = pbuffer.data(idx_dip_hp + 174);

    auto tr_z_yyyyz_y = pbuffer.data(idx_dip_hp + 175);

    auto tr_z_yyyyz_z = pbuffer.data(idx_dip_hp + 176);

#pragma omp simd aligned( \
        pa_y, pa_z, tr_z_yyyy_y, tr_z_yyyyz_x, tr_z_yyyyz_y, tr_z_yyyyz_z, tr_z_yyyz_x, tr_z_yyyz_z, tr_z_yyz_x, tr_z_yyz_z, ts_yyyy_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyyz_x[i] = 3.0 * tr_z_yyz_x[i] * fe_0 + tr_z_yyyz_x[i] * pa_y[i];

        tr_z_yyyyz_y[i] = ts_yyyy_y[i] * fe_0 + tr_z_yyyy_y[i] * pa_z[i];

        tr_z_yyyyz_z[i] = 3.0 * tr_z_yyz_z[i] * fe_0 + tr_z_yyyz_z[i] * pa_y[i];
    }

    // Set up 177-180 components of targeted buffer : HP

    auto tr_z_yyyzz_x = pbuffer.data(idx_dip_hp + 177);

    auto tr_z_yyyzz_y = pbuffer.data(idx_dip_hp + 178);

    auto tr_z_yyyzz_z = pbuffer.data(idx_dip_hp + 179);

#pragma omp simd aligned(pa_y,             \
                             tr_z_yyyzz_x, \
                             tr_z_yyyzz_y, \
                             tr_z_yyyzz_z, \
                             tr_z_yyzz_0,  \
                             tr_z_yyzz_x,  \
                             tr_z_yyzz_y,  \
                             tr_z_yyzz_z,  \
                             tr_z_yzz_x,   \
                             tr_z_yzz_y,   \
                             tr_z_yzz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyzz_x[i] = 2.0 * tr_z_yzz_x[i] * fe_0 + tr_z_yyzz_x[i] * pa_y[i];

        tr_z_yyyzz_y[i] = 2.0 * tr_z_yzz_y[i] * fe_0 + tr_z_yyzz_0[i] * fe_0 + tr_z_yyzz_y[i] * pa_y[i];

        tr_z_yyyzz_z[i] = 2.0 * tr_z_yzz_z[i] * fe_0 + tr_z_yyzz_z[i] * pa_y[i];
    }

    // Set up 180-183 components of targeted buffer : HP

    auto tr_z_yyzzz_x = pbuffer.data(idx_dip_hp + 180);

    auto tr_z_yyzzz_y = pbuffer.data(idx_dip_hp + 181);

    auto tr_z_yyzzz_z = pbuffer.data(idx_dip_hp + 182);

#pragma omp simd aligned(pa_y,             \
                             tr_z_yyzzz_x, \
                             tr_z_yyzzz_y, \
                             tr_z_yyzzz_z, \
                             tr_z_yzzz_0,  \
                             tr_z_yzzz_x,  \
                             tr_z_yzzz_y,  \
                             tr_z_yzzz_z,  \
                             tr_z_zzz_x,   \
                             tr_z_zzz_y,   \
                             tr_z_zzz_z,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyzzz_x[i] = tr_z_zzz_x[i] * fe_0 + tr_z_yzzz_x[i] * pa_y[i];

        tr_z_yyzzz_y[i] = tr_z_zzz_y[i] * fe_0 + tr_z_yzzz_0[i] * fe_0 + tr_z_yzzz_y[i] * pa_y[i];

        tr_z_yyzzz_z[i] = tr_z_zzz_z[i] * fe_0 + tr_z_yzzz_z[i] * pa_y[i];
    }

    // Set up 183-186 components of targeted buffer : HP

    auto tr_z_yzzzz_x = pbuffer.data(idx_dip_hp + 183);

    auto tr_z_yzzzz_y = pbuffer.data(idx_dip_hp + 184);

    auto tr_z_yzzzz_z = pbuffer.data(idx_dip_hp + 185);

#pragma omp simd aligned(pa_y, tr_z_yzzzz_x, tr_z_yzzzz_y, tr_z_yzzzz_z, tr_z_zzzz_0, tr_z_zzzz_x, tr_z_zzzz_y, tr_z_zzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yzzzz_x[i] = tr_z_zzzz_x[i] * pa_y[i];

        tr_z_yzzzz_y[i] = tr_z_zzzz_0[i] * fe_0 + tr_z_zzzz_y[i] * pa_y[i];

        tr_z_yzzzz_z[i] = tr_z_zzzz_z[i] * pa_y[i];
    }

    // Set up 186-189 components of targeted buffer : HP

    auto tr_z_zzzzz_x = pbuffer.data(idx_dip_hp + 186);

    auto tr_z_zzzzz_y = pbuffer.data(idx_dip_hp + 187);

    auto tr_z_zzzzz_z = pbuffer.data(idx_dip_hp + 188);

#pragma omp simd aligned(pa_z,             \
                             tr_z_zzz_x,   \
                             tr_z_zzz_y,   \
                             tr_z_zzz_z,   \
                             tr_z_zzzz_0,  \
                             tr_z_zzzz_x,  \
                             tr_z_zzzz_y,  \
                             tr_z_zzzz_z,  \
                             tr_z_zzzzz_x, \
                             tr_z_zzzzz_y, \
                             tr_z_zzzzz_z, \
                             ts_zzzz_x,    \
                             ts_zzzz_y,    \
                             ts_zzzz_z,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_zzzzz_x[i] = 4.0 * tr_z_zzz_x[i] * fe_0 + ts_zzzz_x[i] * fe_0 + tr_z_zzzz_x[i] * pa_z[i];

        tr_z_zzzzz_y[i] = 4.0 * tr_z_zzz_y[i] * fe_0 + ts_zzzz_y[i] * fe_0 + tr_z_zzzz_y[i] * pa_z[i];

        tr_z_zzzzz_z[i] = 4.0 * tr_z_zzz_z[i] * fe_0 + tr_z_zzzz_0[i] * fe_0 + ts_zzzz_z[i] * fe_0 + tr_z_zzzz_z[i] * pa_z[i];
    }
}

}  // namespace diprec
