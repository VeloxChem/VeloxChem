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

#include "ElectricDipoleMomentumPrimRecDF.hpp"

namespace diprec {  // diprec namespace

auto
comp_prim_electric_dipole_momentum_df(CSimdArray<double>&       pbuffer,
                                      const size_t              idx_dip_df,
                                      const size_t              idx_dip_sf,
                                      const size_t              idx_dip_pd,
                                      const size_t              idx_ovl_pf,
                                      const size_t              idx_dip_pf,
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

    auto tr_x_0_xxx = pbuffer.data(idx_dip_sf);

    auto tr_x_0_xxy = pbuffer.data(idx_dip_sf + 1);

    auto tr_x_0_xxz = pbuffer.data(idx_dip_sf + 2);

    auto tr_x_0_xyy = pbuffer.data(idx_dip_sf + 3);

    auto tr_x_0_xyz = pbuffer.data(idx_dip_sf + 4);

    auto tr_x_0_xzz = pbuffer.data(idx_dip_sf + 5);

    auto tr_x_0_yyy = pbuffer.data(idx_dip_sf + 6);

    auto tr_x_0_yyz = pbuffer.data(idx_dip_sf + 7);

    auto tr_x_0_yzz = pbuffer.data(idx_dip_sf + 8);

    auto tr_x_0_zzz = pbuffer.data(idx_dip_sf + 9);

    auto tr_y_0_xxx = pbuffer.data(idx_dip_sf + 10);

    auto tr_y_0_xxy = pbuffer.data(idx_dip_sf + 11);

    auto tr_y_0_xxz = pbuffer.data(idx_dip_sf + 12);

    auto tr_y_0_xyy = pbuffer.data(idx_dip_sf + 13);

    auto tr_y_0_xyz = pbuffer.data(idx_dip_sf + 14);

    auto tr_y_0_xzz = pbuffer.data(idx_dip_sf + 15);

    auto tr_y_0_yyy = pbuffer.data(idx_dip_sf + 16);

    auto tr_y_0_yyz = pbuffer.data(idx_dip_sf + 17);

    auto tr_y_0_yzz = pbuffer.data(idx_dip_sf + 18);

    auto tr_y_0_zzz = pbuffer.data(idx_dip_sf + 19);

    auto tr_z_0_xxx = pbuffer.data(idx_dip_sf + 20);

    auto tr_z_0_xxy = pbuffer.data(idx_dip_sf + 21);

    auto tr_z_0_xxz = pbuffer.data(idx_dip_sf + 22);

    auto tr_z_0_xyy = pbuffer.data(idx_dip_sf + 23);

    auto tr_z_0_xyz = pbuffer.data(idx_dip_sf + 24);

    auto tr_z_0_xzz = pbuffer.data(idx_dip_sf + 25);

    auto tr_z_0_yyy = pbuffer.data(idx_dip_sf + 26);

    auto tr_z_0_yyz = pbuffer.data(idx_dip_sf + 27);

    auto tr_z_0_yzz = pbuffer.data(idx_dip_sf + 28);

    auto tr_z_0_zzz = pbuffer.data(idx_dip_sf + 29);

    // Set up components of auxiliary buffer : PD

    auto tr_x_x_xx = pbuffer.data(idx_dip_pd);

    auto tr_x_x_xy = pbuffer.data(idx_dip_pd + 1);

    auto tr_x_x_xz = pbuffer.data(idx_dip_pd + 2);

    auto tr_x_x_yy = pbuffer.data(idx_dip_pd + 3);

    auto tr_x_x_yz = pbuffer.data(idx_dip_pd + 4);

    auto tr_x_x_zz = pbuffer.data(idx_dip_pd + 5);

    auto tr_x_y_xx = pbuffer.data(idx_dip_pd + 6);

    auto tr_x_y_xy = pbuffer.data(idx_dip_pd + 7);

    auto tr_x_y_xz = pbuffer.data(idx_dip_pd + 8);

    auto tr_x_y_yy = pbuffer.data(idx_dip_pd + 9);

    auto tr_x_y_yz = pbuffer.data(idx_dip_pd + 10);

    auto tr_x_y_zz = pbuffer.data(idx_dip_pd + 11);

    auto tr_x_z_xx = pbuffer.data(idx_dip_pd + 12);

    auto tr_x_z_xy = pbuffer.data(idx_dip_pd + 13);

    auto tr_x_z_xz = pbuffer.data(idx_dip_pd + 14);

    auto tr_x_z_yy = pbuffer.data(idx_dip_pd + 15);

    auto tr_x_z_yz = pbuffer.data(idx_dip_pd + 16);

    auto tr_x_z_zz = pbuffer.data(idx_dip_pd + 17);

    auto tr_y_x_xx = pbuffer.data(idx_dip_pd + 18);

    auto tr_y_x_xy = pbuffer.data(idx_dip_pd + 19);

    auto tr_y_x_xz = pbuffer.data(idx_dip_pd + 20);

    auto tr_y_x_yy = pbuffer.data(idx_dip_pd + 21);

    auto tr_y_x_yz = pbuffer.data(idx_dip_pd + 22);

    auto tr_y_x_zz = pbuffer.data(idx_dip_pd + 23);

    auto tr_y_y_xx = pbuffer.data(idx_dip_pd + 24);

    auto tr_y_y_xy = pbuffer.data(idx_dip_pd + 25);

    auto tr_y_y_xz = pbuffer.data(idx_dip_pd + 26);

    auto tr_y_y_yy = pbuffer.data(idx_dip_pd + 27);

    auto tr_y_y_yz = pbuffer.data(idx_dip_pd + 28);

    auto tr_y_y_zz = pbuffer.data(idx_dip_pd + 29);

    auto tr_y_z_xx = pbuffer.data(idx_dip_pd + 30);

    auto tr_y_z_xy = pbuffer.data(idx_dip_pd + 31);

    auto tr_y_z_xz = pbuffer.data(idx_dip_pd + 32);

    auto tr_y_z_yy = pbuffer.data(idx_dip_pd + 33);

    auto tr_y_z_yz = pbuffer.data(idx_dip_pd + 34);

    auto tr_y_z_zz = pbuffer.data(idx_dip_pd + 35);

    auto tr_z_x_xx = pbuffer.data(idx_dip_pd + 36);

    auto tr_z_x_xy = pbuffer.data(idx_dip_pd + 37);

    auto tr_z_x_xz = pbuffer.data(idx_dip_pd + 38);

    auto tr_z_x_yy = pbuffer.data(idx_dip_pd + 39);

    auto tr_z_x_yz = pbuffer.data(idx_dip_pd + 40);

    auto tr_z_x_zz = pbuffer.data(idx_dip_pd + 41);

    auto tr_z_y_xx = pbuffer.data(idx_dip_pd + 42);

    auto tr_z_y_xy = pbuffer.data(idx_dip_pd + 43);

    auto tr_z_y_xz = pbuffer.data(idx_dip_pd + 44);

    auto tr_z_y_yy = pbuffer.data(idx_dip_pd + 45);

    auto tr_z_y_yz = pbuffer.data(idx_dip_pd + 46);

    auto tr_z_y_zz = pbuffer.data(idx_dip_pd + 47);

    auto tr_z_z_xx = pbuffer.data(idx_dip_pd + 48);

    auto tr_z_z_xy = pbuffer.data(idx_dip_pd + 49);

    auto tr_z_z_xz = pbuffer.data(idx_dip_pd + 50);

    auto tr_z_z_yy = pbuffer.data(idx_dip_pd + 51);

    auto tr_z_z_yz = pbuffer.data(idx_dip_pd + 52);

    auto tr_z_z_zz = pbuffer.data(idx_dip_pd + 53);

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

    // Set up components of auxiliary buffer : PF

    auto tr_x_x_xxx = pbuffer.data(idx_dip_pf);

    auto tr_x_x_xxy = pbuffer.data(idx_dip_pf + 1);

    auto tr_x_x_xxz = pbuffer.data(idx_dip_pf + 2);

    auto tr_x_x_xyy = pbuffer.data(idx_dip_pf + 3);

    auto tr_x_x_xyz = pbuffer.data(idx_dip_pf + 4);

    auto tr_x_x_xzz = pbuffer.data(idx_dip_pf + 5);

    auto tr_x_x_yyy = pbuffer.data(idx_dip_pf + 6);

    auto tr_x_x_yyz = pbuffer.data(idx_dip_pf + 7);

    auto tr_x_x_yzz = pbuffer.data(idx_dip_pf + 8);

    auto tr_x_x_zzz = pbuffer.data(idx_dip_pf + 9);

    auto tr_x_y_xxx = pbuffer.data(idx_dip_pf + 10);

    auto tr_x_y_xxy = pbuffer.data(idx_dip_pf + 11);

    auto tr_x_y_xxz = pbuffer.data(idx_dip_pf + 12);

    auto tr_x_y_xyy = pbuffer.data(idx_dip_pf + 13);

    auto tr_x_y_xyz = pbuffer.data(idx_dip_pf + 14);

    auto tr_x_y_xzz = pbuffer.data(idx_dip_pf + 15);

    auto tr_x_y_yyy = pbuffer.data(idx_dip_pf + 16);

    auto tr_x_y_yyz = pbuffer.data(idx_dip_pf + 17);

    auto tr_x_y_yzz = pbuffer.data(idx_dip_pf + 18);

    auto tr_x_y_zzz = pbuffer.data(idx_dip_pf + 19);

    auto tr_x_z_xxx = pbuffer.data(idx_dip_pf + 20);

    auto tr_x_z_xxy = pbuffer.data(idx_dip_pf + 21);

    auto tr_x_z_xxz = pbuffer.data(idx_dip_pf + 22);

    auto tr_x_z_xyy = pbuffer.data(idx_dip_pf + 23);

    auto tr_x_z_xyz = pbuffer.data(idx_dip_pf + 24);

    auto tr_x_z_xzz = pbuffer.data(idx_dip_pf + 25);

    auto tr_x_z_yyy = pbuffer.data(idx_dip_pf + 26);

    auto tr_x_z_yyz = pbuffer.data(idx_dip_pf + 27);

    auto tr_x_z_yzz = pbuffer.data(idx_dip_pf + 28);

    auto tr_x_z_zzz = pbuffer.data(idx_dip_pf + 29);

    auto tr_y_x_xxx = pbuffer.data(idx_dip_pf + 30);

    auto tr_y_x_xxy = pbuffer.data(idx_dip_pf + 31);

    auto tr_y_x_xxz = pbuffer.data(idx_dip_pf + 32);

    auto tr_y_x_xyy = pbuffer.data(idx_dip_pf + 33);

    auto tr_y_x_xyz = pbuffer.data(idx_dip_pf + 34);

    auto tr_y_x_xzz = pbuffer.data(idx_dip_pf + 35);

    auto tr_y_x_yyy = pbuffer.data(idx_dip_pf + 36);

    auto tr_y_x_yyz = pbuffer.data(idx_dip_pf + 37);

    auto tr_y_x_yzz = pbuffer.data(idx_dip_pf + 38);

    auto tr_y_x_zzz = pbuffer.data(idx_dip_pf + 39);

    auto tr_y_y_xxx = pbuffer.data(idx_dip_pf + 40);

    auto tr_y_y_xxy = pbuffer.data(idx_dip_pf + 41);

    auto tr_y_y_xxz = pbuffer.data(idx_dip_pf + 42);

    auto tr_y_y_xyy = pbuffer.data(idx_dip_pf + 43);

    auto tr_y_y_xyz = pbuffer.data(idx_dip_pf + 44);

    auto tr_y_y_xzz = pbuffer.data(idx_dip_pf + 45);

    auto tr_y_y_yyy = pbuffer.data(idx_dip_pf + 46);

    auto tr_y_y_yyz = pbuffer.data(idx_dip_pf + 47);

    auto tr_y_y_yzz = pbuffer.data(idx_dip_pf + 48);

    auto tr_y_y_zzz = pbuffer.data(idx_dip_pf + 49);

    auto tr_y_z_xxx = pbuffer.data(idx_dip_pf + 50);

    auto tr_y_z_xxy = pbuffer.data(idx_dip_pf + 51);

    auto tr_y_z_xxz = pbuffer.data(idx_dip_pf + 52);

    auto tr_y_z_xyy = pbuffer.data(idx_dip_pf + 53);

    auto tr_y_z_xyz = pbuffer.data(idx_dip_pf + 54);

    auto tr_y_z_xzz = pbuffer.data(idx_dip_pf + 55);

    auto tr_y_z_yyy = pbuffer.data(idx_dip_pf + 56);

    auto tr_y_z_yyz = pbuffer.data(idx_dip_pf + 57);

    auto tr_y_z_yzz = pbuffer.data(idx_dip_pf + 58);

    auto tr_y_z_zzz = pbuffer.data(idx_dip_pf + 59);

    auto tr_z_x_xxx = pbuffer.data(idx_dip_pf + 60);

    auto tr_z_x_xxy = pbuffer.data(idx_dip_pf + 61);

    auto tr_z_x_xxz = pbuffer.data(idx_dip_pf + 62);

    auto tr_z_x_xyy = pbuffer.data(idx_dip_pf + 63);

    auto tr_z_x_xyz = pbuffer.data(idx_dip_pf + 64);

    auto tr_z_x_xzz = pbuffer.data(idx_dip_pf + 65);

    auto tr_z_x_yyy = pbuffer.data(idx_dip_pf + 66);

    auto tr_z_x_yyz = pbuffer.data(idx_dip_pf + 67);

    auto tr_z_x_yzz = pbuffer.data(idx_dip_pf + 68);

    auto tr_z_x_zzz = pbuffer.data(idx_dip_pf + 69);

    auto tr_z_y_xxx = pbuffer.data(idx_dip_pf + 70);

    auto tr_z_y_xxy = pbuffer.data(idx_dip_pf + 71);

    auto tr_z_y_xxz = pbuffer.data(idx_dip_pf + 72);

    auto tr_z_y_xyy = pbuffer.data(idx_dip_pf + 73);

    auto tr_z_y_xyz = pbuffer.data(idx_dip_pf + 74);

    auto tr_z_y_xzz = pbuffer.data(idx_dip_pf + 75);

    auto tr_z_y_yyy = pbuffer.data(idx_dip_pf + 76);

    auto tr_z_y_yyz = pbuffer.data(idx_dip_pf + 77);

    auto tr_z_y_yzz = pbuffer.data(idx_dip_pf + 78);

    auto tr_z_y_zzz = pbuffer.data(idx_dip_pf + 79);

    auto tr_z_z_xxx = pbuffer.data(idx_dip_pf + 80);

    auto tr_z_z_xxy = pbuffer.data(idx_dip_pf + 81);

    auto tr_z_z_xxz = pbuffer.data(idx_dip_pf + 82);

    auto tr_z_z_xyy = pbuffer.data(idx_dip_pf + 83);

    auto tr_z_z_xyz = pbuffer.data(idx_dip_pf + 84);

    auto tr_z_z_xzz = pbuffer.data(idx_dip_pf + 85);

    auto tr_z_z_yyy = pbuffer.data(idx_dip_pf + 86);

    auto tr_z_z_yyz = pbuffer.data(idx_dip_pf + 87);

    auto tr_z_z_yzz = pbuffer.data(idx_dip_pf + 88);

    auto tr_z_z_zzz = pbuffer.data(idx_dip_pf + 89);

    // Set up 0-10 components of targeted buffer : DF

    auto tr_x_xx_xxx = pbuffer.data(idx_dip_df);

    auto tr_x_xx_xxy = pbuffer.data(idx_dip_df + 1);

    auto tr_x_xx_xxz = pbuffer.data(idx_dip_df + 2);

    auto tr_x_xx_xyy = pbuffer.data(idx_dip_df + 3);

    auto tr_x_xx_xyz = pbuffer.data(idx_dip_df + 4);

    auto tr_x_xx_xzz = pbuffer.data(idx_dip_df + 5);

    auto tr_x_xx_yyy = pbuffer.data(idx_dip_df + 6);

    auto tr_x_xx_yyz = pbuffer.data(idx_dip_df + 7);

    auto tr_x_xx_yzz = pbuffer.data(idx_dip_df + 8);

    auto tr_x_xx_zzz = pbuffer.data(idx_dip_df + 9);

#pragma omp simd aligned(pa_x,            \
                             tr_x_0_xxx,  \
                             tr_x_0_xxy,  \
                             tr_x_0_xxz,  \
                             tr_x_0_xyy,  \
                             tr_x_0_xyz,  \
                             tr_x_0_xzz,  \
                             tr_x_0_yyy,  \
                             tr_x_0_yyz,  \
                             tr_x_0_yzz,  \
                             tr_x_0_zzz,  \
                             tr_x_x_xx,   \
                             tr_x_x_xxx,  \
                             tr_x_x_xxy,  \
                             tr_x_x_xxz,  \
                             tr_x_x_xy,   \
                             tr_x_x_xyy,  \
                             tr_x_x_xyz,  \
                             tr_x_x_xz,   \
                             tr_x_x_xzz,  \
                             tr_x_x_yy,   \
                             tr_x_x_yyy,  \
                             tr_x_x_yyz,  \
                             tr_x_x_yz,   \
                             tr_x_x_yzz,  \
                             tr_x_x_zz,   \
                             tr_x_x_zzz,  \
                             tr_x_xx_xxx, \
                             tr_x_xx_xxy, \
                             tr_x_xx_xxz, \
                             tr_x_xx_xyy, \
                             tr_x_xx_xyz, \
                             tr_x_xx_xzz, \
                             tr_x_xx_yyy, \
                             tr_x_xx_yyz, \
                             tr_x_xx_yzz, \
                             tr_x_xx_zzz, \
                             ts_x_xxx,    \
                             ts_x_xxy,    \
                             ts_x_xxz,    \
                             ts_x_xyy,    \
                             ts_x_xyz,    \
                             ts_x_xzz,    \
                             ts_x_yyy,    \
                             ts_x_yyz,    \
                             ts_x_yzz,    \
                             ts_x_zzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xx_xxx[i] = tr_x_0_xxx[i] * fe_0 + 3.0 * tr_x_x_xx[i] * fe_0 + ts_x_xxx[i] * fe_0 + tr_x_x_xxx[i] * pa_x[i];

        tr_x_xx_xxy[i] = tr_x_0_xxy[i] * fe_0 + 2.0 * tr_x_x_xy[i] * fe_0 + ts_x_xxy[i] * fe_0 + tr_x_x_xxy[i] * pa_x[i];

        tr_x_xx_xxz[i] = tr_x_0_xxz[i] * fe_0 + 2.0 * tr_x_x_xz[i] * fe_0 + ts_x_xxz[i] * fe_0 + tr_x_x_xxz[i] * pa_x[i];

        tr_x_xx_xyy[i] = tr_x_0_xyy[i] * fe_0 + tr_x_x_yy[i] * fe_0 + ts_x_xyy[i] * fe_0 + tr_x_x_xyy[i] * pa_x[i];

        tr_x_xx_xyz[i] = tr_x_0_xyz[i] * fe_0 + tr_x_x_yz[i] * fe_0 + ts_x_xyz[i] * fe_0 + tr_x_x_xyz[i] * pa_x[i];

        tr_x_xx_xzz[i] = tr_x_0_xzz[i] * fe_0 + tr_x_x_zz[i] * fe_0 + ts_x_xzz[i] * fe_0 + tr_x_x_xzz[i] * pa_x[i];

        tr_x_xx_yyy[i] = tr_x_0_yyy[i] * fe_0 + ts_x_yyy[i] * fe_0 + tr_x_x_yyy[i] * pa_x[i];

        tr_x_xx_yyz[i] = tr_x_0_yyz[i] * fe_0 + ts_x_yyz[i] * fe_0 + tr_x_x_yyz[i] * pa_x[i];

        tr_x_xx_yzz[i] = tr_x_0_yzz[i] * fe_0 + ts_x_yzz[i] * fe_0 + tr_x_x_yzz[i] * pa_x[i];

        tr_x_xx_zzz[i] = tr_x_0_zzz[i] * fe_0 + ts_x_zzz[i] * fe_0 + tr_x_x_zzz[i] * pa_x[i];
    }

    // Set up 10-20 components of targeted buffer : DF

    auto tr_x_xy_xxx = pbuffer.data(idx_dip_df + 10);

    auto tr_x_xy_xxy = pbuffer.data(idx_dip_df + 11);

    auto tr_x_xy_xxz = pbuffer.data(idx_dip_df + 12);

    auto tr_x_xy_xyy = pbuffer.data(idx_dip_df + 13);

    auto tr_x_xy_xyz = pbuffer.data(idx_dip_df + 14);

    auto tr_x_xy_xzz = pbuffer.data(idx_dip_df + 15);

    auto tr_x_xy_yyy = pbuffer.data(idx_dip_df + 16);

    auto tr_x_xy_yyz = pbuffer.data(idx_dip_df + 17);

    auto tr_x_xy_yzz = pbuffer.data(idx_dip_df + 18);

    auto tr_x_xy_zzz = pbuffer.data(idx_dip_df + 19);

#pragma omp simd aligned(pa_x,            \
                             pa_y,        \
                             tr_x_x_xx,   \
                             tr_x_x_xxx,  \
                             tr_x_x_xxy,  \
                             tr_x_x_xxz,  \
                             tr_x_x_xy,   \
                             tr_x_x_xyy,  \
                             tr_x_x_xyz,  \
                             tr_x_x_xz,   \
                             tr_x_x_xzz,  \
                             tr_x_x_zzz,  \
                             tr_x_xy_xxx, \
                             tr_x_xy_xxy, \
                             tr_x_xy_xxz, \
                             tr_x_xy_xyy, \
                             tr_x_xy_xyz, \
                             tr_x_xy_xzz, \
                             tr_x_xy_yyy, \
                             tr_x_xy_yyz, \
                             tr_x_xy_yzz, \
                             tr_x_xy_zzz, \
                             tr_x_y_yyy,  \
                             tr_x_y_yyz,  \
                             tr_x_y_yzz,  \
                             ts_y_yyy,    \
                             ts_y_yyz,    \
                             ts_y_yzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xy_xxx[i] = tr_x_x_xxx[i] * pa_y[i];

        tr_x_xy_xxy[i] = tr_x_x_xx[i] * fe_0 + tr_x_x_xxy[i] * pa_y[i];

        tr_x_xy_xxz[i] = tr_x_x_xxz[i] * pa_y[i];

        tr_x_xy_xyy[i] = 2.0 * tr_x_x_xy[i] * fe_0 + tr_x_x_xyy[i] * pa_y[i];

        tr_x_xy_xyz[i] = tr_x_x_xz[i] * fe_0 + tr_x_x_xyz[i] * pa_y[i];

        tr_x_xy_xzz[i] = tr_x_x_xzz[i] * pa_y[i];

        tr_x_xy_yyy[i] = ts_y_yyy[i] * fe_0 + tr_x_y_yyy[i] * pa_x[i];

        tr_x_xy_yyz[i] = ts_y_yyz[i] * fe_0 + tr_x_y_yyz[i] * pa_x[i];

        tr_x_xy_yzz[i] = ts_y_yzz[i] * fe_0 + tr_x_y_yzz[i] * pa_x[i];

        tr_x_xy_zzz[i] = tr_x_x_zzz[i] * pa_y[i];
    }

    // Set up 20-30 components of targeted buffer : DF

    auto tr_x_xz_xxx = pbuffer.data(idx_dip_df + 20);

    auto tr_x_xz_xxy = pbuffer.data(idx_dip_df + 21);

    auto tr_x_xz_xxz = pbuffer.data(idx_dip_df + 22);

    auto tr_x_xz_xyy = pbuffer.data(idx_dip_df + 23);

    auto tr_x_xz_xyz = pbuffer.data(idx_dip_df + 24);

    auto tr_x_xz_xzz = pbuffer.data(idx_dip_df + 25);

    auto tr_x_xz_yyy = pbuffer.data(idx_dip_df + 26);

    auto tr_x_xz_yyz = pbuffer.data(idx_dip_df + 27);

    auto tr_x_xz_yzz = pbuffer.data(idx_dip_df + 28);

    auto tr_x_xz_zzz = pbuffer.data(idx_dip_df + 29);

#pragma omp simd aligned(pa_x,            \
                             pa_z,        \
                             tr_x_x_xx,   \
                             tr_x_x_xxx,  \
                             tr_x_x_xxy,  \
                             tr_x_x_xxz,  \
                             tr_x_x_xy,   \
                             tr_x_x_xyy,  \
                             tr_x_x_xyz,  \
                             tr_x_x_xz,   \
                             tr_x_x_xzz,  \
                             tr_x_x_yyy,  \
                             tr_x_xz_xxx, \
                             tr_x_xz_xxy, \
                             tr_x_xz_xxz, \
                             tr_x_xz_xyy, \
                             tr_x_xz_xyz, \
                             tr_x_xz_xzz, \
                             tr_x_xz_yyy, \
                             tr_x_xz_yyz, \
                             tr_x_xz_yzz, \
                             tr_x_xz_zzz, \
                             tr_x_z_yyz,  \
                             tr_x_z_yzz,  \
                             tr_x_z_zzz,  \
                             ts_z_yyz,    \
                             ts_z_yzz,    \
                             ts_z_zzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xz_xxx[i] = tr_x_x_xxx[i] * pa_z[i];

        tr_x_xz_xxy[i] = tr_x_x_xxy[i] * pa_z[i];

        tr_x_xz_xxz[i] = tr_x_x_xx[i] * fe_0 + tr_x_x_xxz[i] * pa_z[i];

        tr_x_xz_xyy[i] = tr_x_x_xyy[i] * pa_z[i];

        tr_x_xz_xyz[i] = tr_x_x_xy[i] * fe_0 + tr_x_x_xyz[i] * pa_z[i];

        tr_x_xz_xzz[i] = 2.0 * tr_x_x_xz[i] * fe_0 + tr_x_x_xzz[i] * pa_z[i];

        tr_x_xz_yyy[i] = tr_x_x_yyy[i] * pa_z[i];

        tr_x_xz_yyz[i] = ts_z_yyz[i] * fe_0 + tr_x_z_yyz[i] * pa_x[i];

        tr_x_xz_yzz[i] = ts_z_yzz[i] * fe_0 + tr_x_z_yzz[i] * pa_x[i];

        tr_x_xz_zzz[i] = ts_z_zzz[i] * fe_0 + tr_x_z_zzz[i] * pa_x[i];
    }

    // Set up 30-40 components of targeted buffer : DF

    auto tr_x_yy_xxx = pbuffer.data(idx_dip_df + 30);

    auto tr_x_yy_xxy = pbuffer.data(idx_dip_df + 31);

    auto tr_x_yy_xxz = pbuffer.data(idx_dip_df + 32);

    auto tr_x_yy_xyy = pbuffer.data(idx_dip_df + 33);

    auto tr_x_yy_xyz = pbuffer.data(idx_dip_df + 34);

    auto tr_x_yy_xzz = pbuffer.data(idx_dip_df + 35);

    auto tr_x_yy_yyy = pbuffer.data(idx_dip_df + 36);

    auto tr_x_yy_yyz = pbuffer.data(idx_dip_df + 37);

    auto tr_x_yy_yzz = pbuffer.data(idx_dip_df + 38);

    auto tr_x_yy_zzz = pbuffer.data(idx_dip_df + 39);

#pragma omp simd aligned(pa_y,            \
                             tr_x_0_xxx,  \
                             tr_x_0_xxy,  \
                             tr_x_0_xxz,  \
                             tr_x_0_xyy,  \
                             tr_x_0_xyz,  \
                             tr_x_0_xzz,  \
                             tr_x_0_yyy,  \
                             tr_x_0_yyz,  \
                             tr_x_0_yzz,  \
                             tr_x_0_zzz,  \
                             tr_x_y_xx,   \
                             tr_x_y_xxx,  \
                             tr_x_y_xxy,  \
                             tr_x_y_xxz,  \
                             tr_x_y_xy,   \
                             tr_x_y_xyy,  \
                             tr_x_y_xyz,  \
                             tr_x_y_xz,   \
                             tr_x_y_xzz,  \
                             tr_x_y_yy,   \
                             tr_x_y_yyy,  \
                             tr_x_y_yyz,  \
                             tr_x_y_yz,   \
                             tr_x_y_yzz,  \
                             tr_x_y_zz,   \
                             tr_x_y_zzz,  \
                             tr_x_yy_xxx, \
                             tr_x_yy_xxy, \
                             tr_x_yy_xxz, \
                             tr_x_yy_xyy, \
                             tr_x_yy_xyz, \
                             tr_x_yy_xzz, \
                             tr_x_yy_yyy, \
                             tr_x_yy_yyz, \
                             tr_x_yy_yzz, \
                             tr_x_yy_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yy_xxx[i] = tr_x_0_xxx[i] * fe_0 + tr_x_y_xxx[i] * pa_y[i];

        tr_x_yy_xxy[i] = tr_x_0_xxy[i] * fe_0 + tr_x_y_xx[i] * fe_0 + tr_x_y_xxy[i] * pa_y[i];

        tr_x_yy_xxz[i] = tr_x_0_xxz[i] * fe_0 + tr_x_y_xxz[i] * pa_y[i];

        tr_x_yy_xyy[i] = tr_x_0_xyy[i] * fe_0 + 2.0 * tr_x_y_xy[i] * fe_0 + tr_x_y_xyy[i] * pa_y[i];

        tr_x_yy_xyz[i] = tr_x_0_xyz[i] * fe_0 + tr_x_y_xz[i] * fe_0 + tr_x_y_xyz[i] * pa_y[i];

        tr_x_yy_xzz[i] = tr_x_0_xzz[i] * fe_0 + tr_x_y_xzz[i] * pa_y[i];

        tr_x_yy_yyy[i] = tr_x_0_yyy[i] * fe_0 + 3.0 * tr_x_y_yy[i] * fe_0 + tr_x_y_yyy[i] * pa_y[i];

        tr_x_yy_yyz[i] = tr_x_0_yyz[i] * fe_0 + 2.0 * tr_x_y_yz[i] * fe_0 + tr_x_y_yyz[i] * pa_y[i];

        tr_x_yy_yzz[i] = tr_x_0_yzz[i] * fe_0 + tr_x_y_zz[i] * fe_0 + tr_x_y_yzz[i] * pa_y[i];

        tr_x_yy_zzz[i] = tr_x_0_zzz[i] * fe_0 + tr_x_y_zzz[i] * pa_y[i];
    }

    // Set up 40-50 components of targeted buffer : DF

    auto tr_x_yz_xxx = pbuffer.data(idx_dip_df + 40);

    auto tr_x_yz_xxy = pbuffer.data(idx_dip_df + 41);

    auto tr_x_yz_xxz = pbuffer.data(idx_dip_df + 42);

    auto tr_x_yz_xyy = pbuffer.data(idx_dip_df + 43);

    auto tr_x_yz_xyz = pbuffer.data(idx_dip_df + 44);

    auto tr_x_yz_xzz = pbuffer.data(idx_dip_df + 45);

    auto tr_x_yz_yyy = pbuffer.data(idx_dip_df + 46);

    auto tr_x_yz_yyz = pbuffer.data(idx_dip_df + 47);

    auto tr_x_yz_yzz = pbuffer.data(idx_dip_df + 48);

    auto tr_x_yz_zzz = pbuffer.data(idx_dip_df + 49);

#pragma omp simd aligned(pa_y,            \
                             pa_z,        \
                             tr_x_y_xxy,  \
                             tr_x_y_xyy,  \
                             tr_x_y_yyy,  \
                             tr_x_yz_xxx, \
                             tr_x_yz_xxy, \
                             tr_x_yz_xxz, \
                             tr_x_yz_xyy, \
                             tr_x_yz_xyz, \
                             tr_x_yz_xzz, \
                             tr_x_yz_yyy, \
                             tr_x_yz_yyz, \
                             tr_x_yz_yzz, \
                             tr_x_yz_zzz, \
                             tr_x_z_xxx,  \
                             tr_x_z_xxz,  \
                             tr_x_z_xyz,  \
                             tr_x_z_xz,   \
                             tr_x_z_xzz,  \
                             tr_x_z_yyz,  \
                             tr_x_z_yz,   \
                             tr_x_z_yzz,  \
                             tr_x_z_zz,   \
                             tr_x_z_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yz_xxx[i] = tr_x_z_xxx[i] * pa_y[i];

        tr_x_yz_xxy[i] = tr_x_y_xxy[i] * pa_z[i];

        tr_x_yz_xxz[i] = tr_x_z_xxz[i] * pa_y[i];

        tr_x_yz_xyy[i] = tr_x_y_xyy[i] * pa_z[i];

        tr_x_yz_xyz[i] = tr_x_z_xz[i] * fe_0 + tr_x_z_xyz[i] * pa_y[i];

        tr_x_yz_xzz[i] = tr_x_z_xzz[i] * pa_y[i];

        tr_x_yz_yyy[i] = tr_x_y_yyy[i] * pa_z[i];

        tr_x_yz_yyz[i] = 2.0 * tr_x_z_yz[i] * fe_0 + tr_x_z_yyz[i] * pa_y[i];

        tr_x_yz_yzz[i] = tr_x_z_zz[i] * fe_0 + tr_x_z_yzz[i] * pa_y[i];

        tr_x_yz_zzz[i] = tr_x_z_zzz[i] * pa_y[i];
    }

    // Set up 50-60 components of targeted buffer : DF

    auto tr_x_zz_xxx = pbuffer.data(idx_dip_df + 50);

    auto tr_x_zz_xxy = pbuffer.data(idx_dip_df + 51);

    auto tr_x_zz_xxz = pbuffer.data(idx_dip_df + 52);

    auto tr_x_zz_xyy = pbuffer.data(idx_dip_df + 53);

    auto tr_x_zz_xyz = pbuffer.data(idx_dip_df + 54);

    auto tr_x_zz_xzz = pbuffer.data(idx_dip_df + 55);

    auto tr_x_zz_yyy = pbuffer.data(idx_dip_df + 56);

    auto tr_x_zz_yyz = pbuffer.data(idx_dip_df + 57);

    auto tr_x_zz_yzz = pbuffer.data(idx_dip_df + 58);

    auto tr_x_zz_zzz = pbuffer.data(idx_dip_df + 59);

#pragma omp simd aligned(pa_z,            \
                             tr_x_0_xxx,  \
                             tr_x_0_xxy,  \
                             tr_x_0_xxz,  \
                             tr_x_0_xyy,  \
                             tr_x_0_xyz,  \
                             tr_x_0_xzz,  \
                             tr_x_0_yyy,  \
                             tr_x_0_yyz,  \
                             tr_x_0_yzz,  \
                             tr_x_0_zzz,  \
                             tr_x_z_xx,   \
                             tr_x_z_xxx,  \
                             tr_x_z_xxy,  \
                             tr_x_z_xxz,  \
                             tr_x_z_xy,   \
                             tr_x_z_xyy,  \
                             tr_x_z_xyz,  \
                             tr_x_z_xz,   \
                             tr_x_z_xzz,  \
                             tr_x_z_yy,   \
                             tr_x_z_yyy,  \
                             tr_x_z_yyz,  \
                             tr_x_z_yz,   \
                             tr_x_z_yzz,  \
                             tr_x_z_zz,   \
                             tr_x_z_zzz,  \
                             tr_x_zz_xxx, \
                             tr_x_zz_xxy, \
                             tr_x_zz_xxz, \
                             tr_x_zz_xyy, \
                             tr_x_zz_xyz, \
                             tr_x_zz_xzz, \
                             tr_x_zz_yyy, \
                             tr_x_zz_yyz, \
                             tr_x_zz_yzz, \
                             tr_x_zz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_zz_xxx[i] = tr_x_0_xxx[i] * fe_0 + tr_x_z_xxx[i] * pa_z[i];

        tr_x_zz_xxy[i] = tr_x_0_xxy[i] * fe_0 + tr_x_z_xxy[i] * pa_z[i];

        tr_x_zz_xxz[i] = tr_x_0_xxz[i] * fe_0 + tr_x_z_xx[i] * fe_0 + tr_x_z_xxz[i] * pa_z[i];

        tr_x_zz_xyy[i] = tr_x_0_xyy[i] * fe_0 + tr_x_z_xyy[i] * pa_z[i];

        tr_x_zz_xyz[i] = tr_x_0_xyz[i] * fe_0 + tr_x_z_xy[i] * fe_0 + tr_x_z_xyz[i] * pa_z[i];

        tr_x_zz_xzz[i] = tr_x_0_xzz[i] * fe_0 + 2.0 * tr_x_z_xz[i] * fe_0 + tr_x_z_xzz[i] * pa_z[i];

        tr_x_zz_yyy[i] = tr_x_0_yyy[i] * fe_0 + tr_x_z_yyy[i] * pa_z[i];

        tr_x_zz_yyz[i] = tr_x_0_yyz[i] * fe_0 + tr_x_z_yy[i] * fe_0 + tr_x_z_yyz[i] * pa_z[i];

        tr_x_zz_yzz[i] = tr_x_0_yzz[i] * fe_0 + 2.0 * tr_x_z_yz[i] * fe_0 + tr_x_z_yzz[i] * pa_z[i];

        tr_x_zz_zzz[i] = tr_x_0_zzz[i] * fe_0 + 3.0 * tr_x_z_zz[i] * fe_0 + tr_x_z_zzz[i] * pa_z[i];
    }

    // Set up 60-70 components of targeted buffer : DF

    auto tr_y_xx_xxx = pbuffer.data(idx_dip_df + 60);

    auto tr_y_xx_xxy = pbuffer.data(idx_dip_df + 61);

    auto tr_y_xx_xxz = pbuffer.data(idx_dip_df + 62);

    auto tr_y_xx_xyy = pbuffer.data(idx_dip_df + 63);

    auto tr_y_xx_xyz = pbuffer.data(idx_dip_df + 64);

    auto tr_y_xx_xzz = pbuffer.data(idx_dip_df + 65);

    auto tr_y_xx_yyy = pbuffer.data(idx_dip_df + 66);

    auto tr_y_xx_yyz = pbuffer.data(idx_dip_df + 67);

    auto tr_y_xx_yzz = pbuffer.data(idx_dip_df + 68);

    auto tr_y_xx_zzz = pbuffer.data(idx_dip_df + 69);

#pragma omp simd aligned(pa_x,            \
                             tr_y_0_xxx,  \
                             tr_y_0_xxy,  \
                             tr_y_0_xxz,  \
                             tr_y_0_xyy,  \
                             tr_y_0_xyz,  \
                             tr_y_0_xzz,  \
                             tr_y_0_yyy,  \
                             tr_y_0_yyz,  \
                             tr_y_0_yzz,  \
                             tr_y_0_zzz,  \
                             tr_y_x_xx,   \
                             tr_y_x_xxx,  \
                             tr_y_x_xxy,  \
                             tr_y_x_xxz,  \
                             tr_y_x_xy,   \
                             tr_y_x_xyy,  \
                             tr_y_x_xyz,  \
                             tr_y_x_xz,   \
                             tr_y_x_xzz,  \
                             tr_y_x_yy,   \
                             tr_y_x_yyy,  \
                             tr_y_x_yyz,  \
                             tr_y_x_yz,   \
                             tr_y_x_yzz,  \
                             tr_y_x_zz,   \
                             tr_y_x_zzz,  \
                             tr_y_xx_xxx, \
                             tr_y_xx_xxy, \
                             tr_y_xx_xxz, \
                             tr_y_xx_xyy, \
                             tr_y_xx_xyz, \
                             tr_y_xx_xzz, \
                             tr_y_xx_yyy, \
                             tr_y_xx_yyz, \
                             tr_y_xx_yzz, \
                             tr_y_xx_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xx_xxx[i] = tr_y_0_xxx[i] * fe_0 + 3.0 * tr_y_x_xx[i] * fe_0 + tr_y_x_xxx[i] * pa_x[i];

        tr_y_xx_xxy[i] = tr_y_0_xxy[i] * fe_0 + 2.0 * tr_y_x_xy[i] * fe_0 + tr_y_x_xxy[i] * pa_x[i];

        tr_y_xx_xxz[i] = tr_y_0_xxz[i] * fe_0 + 2.0 * tr_y_x_xz[i] * fe_0 + tr_y_x_xxz[i] * pa_x[i];

        tr_y_xx_xyy[i] = tr_y_0_xyy[i] * fe_0 + tr_y_x_yy[i] * fe_0 + tr_y_x_xyy[i] * pa_x[i];

        tr_y_xx_xyz[i] = tr_y_0_xyz[i] * fe_0 + tr_y_x_yz[i] * fe_0 + tr_y_x_xyz[i] * pa_x[i];

        tr_y_xx_xzz[i] = tr_y_0_xzz[i] * fe_0 + tr_y_x_zz[i] * fe_0 + tr_y_x_xzz[i] * pa_x[i];

        tr_y_xx_yyy[i] = tr_y_0_yyy[i] * fe_0 + tr_y_x_yyy[i] * pa_x[i];

        tr_y_xx_yyz[i] = tr_y_0_yyz[i] * fe_0 + tr_y_x_yyz[i] * pa_x[i];

        tr_y_xx_yzz[i] = tr_y_0_yzz[i] * fe_0 + tr_y_x_yzz[i] * pa_x[i];

        tr_y_xx_zzz[i] = tr_y_0_zzz[i] * fe_0 + tr_y_x_zzz[i] * pa_x[i];
    }

    // Set up 70-80 components of targeted buffer : DF

    auto tr_y_xy_xxx = pbuffer.data(idx_dip_df + 70);

    auto tr_y_xy_xxy = pbuffer.data(idx_dip_df + 71);

    auto tr_y_xy_xxz = pbuffer.data(idx_dip_df + 72);

    auto tr_y_xy_xyy = pbuffer.data(idx_dip_df + 73);

    auto tr_y_xy_xyz = pbuffer.data(idx_dip_df + 74);

    auto tr_y_xy_xzz = pbuffer.data(idx_dip_df + 75);

    auto tr_y_xy_yyy = pbuffer.data(idx_dip_df + 76);

    auto tr_y_xy_yyz = pbuffer.data(idx_dip_df + 77);

    auto tr_y_xy_yzz = pbuffer.data(idx_dip_df + 78);

    auto tr_y_xy_zzz = pbuffer.data(idx_dip_df + 79);

#pragma omp simd aligned(pa_x,            \
                             tr_y_xy_xxx, \
                             tr_y_xy_xxy, \
                             tr_y_xy_xxz, \
                             tr_y_xy_xyy, \
                             tr_y_xy_xyz, \
                             tr_y_xy_xzz, \
                             tr_y_xy_yyy, \
                             tr_y_xy_yyz, \
                             tr_y_xy_yzz, \
                             tr_y_xy_zzz, \
                             tr_y_y_xx,   \
                             tr_y_y_xxx,  \
                             tr_y_y_xxy,  \
                             tr_y_y_xxz,  \
                             tr_y_y_xy,   \
                             tr_y_y_xyy,  \
                             tr_y_y_xyz,  \
                             tr_y_y_xz,   \
                             tr_y_y_xzz,  \
                             tr_y_y_yy,   \
                             tr_y_y_yyy,  \
                             tr_y_y_yyz,  \
                             tr_y_y_yz,   \
                             tr_y_y_yzz,  \
                             tr_y_y_zz,   \
                             tr_y_y_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xy_xxx[i] = 3.0 * tr_y_y_xx[i] * fe_0 + tr_y_y_xxx[i] * pa_x[i];

        tr_y_xy_xxy[i] = 2.0 * tr_y_y_xy[i] * fe_0 + tr_y_y_xxy[i] * pa_x[i];

        tr_y_xy_xxz[i] = 2.0 * tr_y_y_xz[i] * fe_0 + tr_y_y_xxz[i] * pa_x[i];

        tr_y_xy_xyy[i] = tr_y_y_yy[i] * fe_0 + tr_y_y_xyy[i] * pa_x[i];

        tr_y_xy_xyz[i] = tr_y_y_yz[i] * fe_0 + tr_y_y_xyz[i] * pa_x[i];

        tr_y_xy_xzz[i] = tr_y_y_zz[i] * fe_0 + tr_y_y_xzz[i] * pa_x[i];

        tr_y_xy_yyy[i] = tr_y_y_yyy[i] * pa_x[i];

        tr_y_xy_yyz[i] = tr_y_y_yyz[i] * pa_x[i];

        tr_y_xy_yzz[i] = tr_y_y_yzz[i] * pa_x[i];

        tr_y_xy_zzz[i] = tr_y_y_zzz[i] * pa_x[i];
    }

    // Set up 80-90 components of targeted buffer : DF

    auto tr_y_xz_xxx = pbuffer.data(idx_dip_df + 80);

    auto tr_y_xz_xxy = pbuffer.data(idx_dip_df + 81);

    auto tr_y_xz_xxz = pbuffer.data(idx_dip_df + 82);

    auto tr_y_xz_xyy = pbuffer.data(idx_dip_df + 83);

    auto tr_y_xz_xyz = pbuffer.data(idx_dip_df + 84);

    auto tr_y_xz_xzz = pbuffer.data(idx_dip_df + 85);

    auto tr_y_xz_yyy = pbuffer.data(idx_dip_df + 86);

    auto tr_y_xz_yyz = pbuffer.data(idx_dip_df + 87);

    auto tr_y_xz_yzz = pbuffer.data(idx_dip_df + 88);

    auto tr_y_xz_zzz = pbuffer.data(idx_dip_df + 89);

#pragma omp simd aligned(pa_x,            \
                             pa_z,        \
                             tr_y_x_xxx,  \
                             tr_y_x_xxy,  \
                             tr_y_x_xyy,  \
                             tr_y_xz_xxx, \
                             tr_y_xz_xxy, \
                             tr_y_xz_xxz, \
                             tr_y_xz_xyy, \
                             tr_y_xz_xyz, \
                             tr_y_xz_xzz, \
                             tr_y_xz_yyy, \
                             tr_y_xz_yyz, \
                             tr_y_xz_yzz, \
                             tr_y_xz_zzz, \
                             tr_y_z_xxz,  \
                             tr_y_z_xyz,  \
                             tr_y_z_xz,   \
                             tr_y_z_xzz,  \
                             tr_y_z_yyy,  \
                             tr_y_z_yyz,  \
                             tr_y_z_yz,   \
                             tr_y_z_yzz,  \
                             tr_y_z_zz,   \
                             tr_y_z_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xz_xxx[i] = tr_y_x_xxx[i] * pa_z[i];

        tr_y_xz_xxy[i] = tr_y_x_xxy[i] * pa_z[i];

        tr_y_xz_xxz[i] = 2.0 * tr_y_z_xz[i] * fe_0 + tr_y_z_xxz[i] * pa_x[i];

        tr_y_xz_xyy[i] = tr_y_x_xyy[i] * pa_z[i];

        tr_y_xz_xyz[i] = tr_y_z_yz[i] * fe_0 + tr_y_z_xyz[i] * pa_x[i];

        tr_y_xz_xzz[i] = tr_y_z_zz[i] * fe_0 + tr_y_z_xzz[i] * pa_x[i];

        tr_y_xz_yyy[i] = tr_y_z_yyy[i] * pa_x[i];

        tr_y_xz_yyz[i] = tr_y_z_yyz[i] * pa_x[i];

        tr_y_xz_yzz[i] = tr_y_z_yzz[i] * pa_x[i];

        tr_y_xz_zzz[i] = tr_y_z_zzz[i] * pa_x[i];
    }

    // Set up 90-100 components of targeted buffer : DF

    auto tr_y_yy_xxx = pbuffer.data(idx_dip_df + 90);

    auto tr_y_yy_xxy = pbuffer.data(idx_dip_df + 91);

    auto tr_y_yy_xxz = pbuffer.data(idx_dip_df + 92);

    auto tr_y_yy_xyy = pbuffer.data(idx_dip_df + 93);

    auto tr_y_yy_xyz = pbuffer.data(idx_dip_df + 94);

    auto tr_y_yy_xzz = pbuffer.data(idx_dip_df + 95);

    auto tr_y_yy_yyy = pbuffer.data(idx_dip_df + 96);

    auto tr_y_yy_yyz = pbuffer.data(idx_dip_df + 97);

    auto tr_y_yy_yzz = pbuffer.data(idx_dip_df + 98);

    auto tr_y_yy_zzz = pbuffer.data(idx_dip_df + 99);

#pragma omp simd aligned(pa_y,            \
                             tr_y_0_xxx,  \
                             tr_y_0_xxy,  \
                             tr_y_0_xxz,  \
                             tr_y_0_xyy,  \
                             tr_y_0_xyz,  \
                             tr_y_0_xzz,  \
                             tr_y_0_yyy,  \
                             tr_y_0_yyz,  \
                             tr_y_0_yzz,  \
                             tr_y_0_zzz,  \
                             tr_y_y_xx,   \
                             tr_y_y_xxx,  \
                             tr_y_y_xxy,  \
                             tr_y_y_xxz,  \
                             tr_y_y_xy,   \
                             tr_y_y_xyy,  \
                             tr_y_y_xyz,  \
                             tr_y_y_xz,   \
                             tr_y_y_xzz,  \
                             tr_y_y_yy,   \
                             tr_y_y_yyy,  \
                             tr_y_y_yyz,  \
                             tr_y_y_yz,   \
                             tr_y_y_yzz,  \
                             tr_y_y_zz,   \
                             tr_y_y_zzz,  \
                             tr_y_yy_xxx, \
                             tr_y_yy_xxy, \
                             tr_y_yy_xxz, \
                             tr_y_yy_xyy, \
                             tr_y_yy_xyz, \
                             tr_y_yy_xzz, \
                             tr_y_yy_yyy, \
                             tr_y_yy_yyz, \
                             tr_y_yy_yzz, \
                             tr_y_yy_zzz, \
                             ts_y_xxx,    \
                             ts_y_xxy,    \
                             ts_y_xxz,    \
                             ts_y_xyy,    \
                             ts_y_xyz,    \
                             ts_y_xzz,    \
                             ts_y_yyy,    \
                             ts_y_yyz,    \
                             ts_y_yzz,    \
                             ts_y_zzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yy_xxx[i] = tr_y_0_xxx[i] * fe_0 + ts_y_xxx[i] * fe_0 + tr_y_y_xxx[i] * pa_y[i];

        tr_y_yy_xxy[i] = tr_y_0_xxy[i] * fe_0 + tr_y_y_xx[i] * fe_0 + ts_y_xxy[i] * fe_0 + tr_y_y_xxy[i] * pa_y[i];

        tr_y_yy_xxz[i] = tr_y_0_xxz[i] * fe_0 + ts_y_xxz[i] * fe_0 + tr_y_y_xxz[i] * pa_y[i];

        tr_y_yy_xyy[i] = tr_y_0_xyy[i] * fe_0 + 2.0 * tr_y_y_xy[i] * fe_0 + ts_y_xyy[i] * fe_0 + tr_y_y_xyy[i] * pa_y[i];

        tr_y_yy_xyz[i] = tr_y_0_xyz[i] * fe_0 + tr_y_y_xz[i] * fe_0 + ts_y_xyz[i] * fe_0 + tr_y_y_xyz[i] * pa_y[i];

        tr_y_yy_xzz[i] = tr_y_0_xzz[i] * fe_0 + ts_y_xzz[i] * fe_0 + tr_y_y_xzz[i] * pa_y[i];

        tr_y_yy_yyy[i] = tr_y_0_yyy[i] * fe_0 + 3.0 * tr_y_y_yy[i] * fe_0 + ts_y_yyy[i] * fe_0 + tr_y_y_yyy[i] * pa_y[i];

        tr_y_yy_yyz[i] = tr_y_0_yyz[i] * fe_0 + 2.0 * tr_y_y_yz[i] * fe_0 + ts_y_yyz[i] * fe_0 + tr_y_y_yyz[i] * pa_y[i];

        tr_y_yy_yzz[i] = tr_y_0_yzz[i] * fe_0 + tr_y_y_zz[i] * fe_0 + ts_y_yzz[i] * fe_0 + tr_y_y_yzz[i] * pa_y[i];

        tr_y_yy_zzz[i] = tr_y_0_zzz[i] * fe_0 + ts_y_zzz[i] * fe_0 + tr_y_y_zzz[i] * pa_y[i];
    }

    // Set up 100-110 components of targeted buffer : DF

    auto tr_y_yz_xxx = pbuffer.data(idx_dip_df + 100);

    auto tr_y_yz_xxy = pbuffer.data(idx_dip_df + 101);

    auto tr_y_yz_xxz = pbuffer.data(idx_dip_df + 102);

    auto tr_y_yz_xyy = pbuffer.data(idx_dip_df + 103);

    auto tr_y_yz_xyz = pbuffer.data(idx_dip_df + 104);

    auto tr_y_yz_xzz = pbuffer.data(idx_dip_df + 105);

    auto tr_y_yz_yyy = pbuffer.data(idx_dip_df + 106);

    auto tr_y_yz_yyz = pbuffer.data(idx_dip_df + 107);

    auto tr_y_yz_yzz = pbuffer.data(idx_dip_df + 108);

    auto tr_y_yz_zzz = pbuffer.data(idx_dip_df + 109);

#pragma omp simd aligned(pa_y,            \
                             pa_z,        \
                             tr_y_y_xxx,  \
                             tr_y_y_xxy,  \
                             tr_y_y_xy,   \
                             tr_y_y_xyy,  \
                             tr_y_y_xyz,  \
                             tr_y_y_yy,   \
                             tr_y_y_yyy,  \
                             tr_y_y_yyz,  \
                             tr_y_y_yz,   \
                             tr_y_y_yzz,  \
                             tr_y_yz_xxx, \
                             tr_y_yz_xxy, \
                             tr_y_yz_xxz, \
                             tr_y_yz_xyy, \
                             tr_y_yz_xyz, \
                             tr_y_yz_xzz, \
                             tr_y_yz_yyy, \
                             tr_y_yz_yyz, \
                             tr_y_yz_yzz, \
                             tr_y_yz_zzz, \
                             tr_y_z_xxz,  \
                             tr_y_z_xzz,  \
                             tr_y_z_zzz,  \
                             ts_z_xxz,    \
                             ts_z_xzz,    \
                             ts_z_zzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yz_xxx[i] = tr_y_y_xxx[i] * pa_z[i];

        tr_y_yz_xxy[i] = tr_y_y_xxy[i] * pa_z[i];

        tr_y_yz_xxz[i] = ts_z_xxz[i] * fe_0 + tr_y_z_xxz[i] * pa_y[i];

        tr_y_yz_xyy[i] = tr_y_y_xyy[i] * pa_z[i];

        tr_y_yz_xyz[i] = tr_y_y_xy[i] * fe_0 + tr_y_y_xyz[i] * pa_z[i];

        tr_y_yz_xzz[i] = ts_z_xzz[i] * fe_0 + tr_y_z_xzz[i] * pa_y[i];

        tr_y_yz_yyy[i] = tr_y_y_yyy[i] * pa_z[i];

        tr_y_yz_yyz[i] = tr_y_y_yy[i] * fe_0 + tr_y_y_yyz[i] * pa_z[i];

        tr_y_yz_yzz[i] = 2.0 * tr_y_y_yz[i] * fe_0 + tr_y_y_yzz[i] * pa_z[i];

        tr_y_yz_zzz[i] = ts_z_zzz[i] * fe_0 + tr_y_z_zzz[i] * pa_y[i];
    }

    // Set up 110-120 components of targeted buffer : DF

    auto tr_y_zz_xxx = pbuffer.data(idx_dip_df + 110);

    auto tr_y_zz_xxy = pbuffer.data(idx_dip_df + 111);

    auto tr_y_zz_xxz = pbuffer.data(idx_dip_df + 112);

    auto tr_y_zz_xyy = pbuffer.data(idx_dip_df + 113);

    auto tr_y_zz_xyz = pbuffer.data(idx_dip_df + 114);

    auto tr_y_zz_xzz = pbuffer.data(idx_dip_df + 115);

    auto tr_y_zz_yyy = pbuffer.data(idx_dip_df + 116);

    auto tr_y_zz_yyz = pbuffer.data(idx_dip_df + 117);

    auto tr_y_zz_yzz = pbuffer.data(idx_dip_df + 118);

    auto tr_y_zz_zzz = pbuffer.data(idx_dip_df + 119);

#pragma omp simd aligned(pa_z,            \
                             tr_y_0_xxx,  \
                             tr_y_0_xxy,  \
                             tr_y_0_xxz,  \
                             tr_y_0_xyy,  \
                             tr_y_0_xyz,  \
                             tr_y_0_xzz,  \
                             tr_y_0_yyy,  \
                             tr_y_0_yyz,  \
                             tr_y_0_yzz,  \
                             tr_y_0_zzz,  \
                             tr_y_z_xx,   \
                             tr_y_z_xxx,  \
                             tr_y_z_xxy,  \
                             tr_y_z_xxz,  \
                             tr_y_z_xy,   \
                             tr_y_z_xyy,  \
                             tr_y_z_xyz,  \
                             tr_y_z_xz,   \
                             tr_y_z_xzz,  \
                             tr_y_z_yy,   \
                             tr_y_z_yyy,  \
                             tr_y_z_yyz,  \
                             tr_y_z_yz,   \
                             tr_y_z_yzz,  \
                             tr_y_z_zz,   \
                             tr_y_z_zzz,  \
                             tr_y_zz_xxx, \
                             tr_y_zz_xxy, \
                             tr_y_zz_xxz, \
                             tr_y_zz_xyy, \
                             tr_y_zz_xyz, \
                             tr_y_zz_xzz, \
                             tr_y_zz_yyy, \
                             tr_y_zz_yyz, \
                             tr_y_zz_yzz, \
                             tr_y_zz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_zz_xxx[i] = tr_y_0_xxx[i] * fe_0 + tr_y_z_xxx[i] * pa_z[i];

        tr_y_zz_xxy[i] = tr_y_0_xxy[i] * fe_0 + tr_y_z_xxy[i] * pa_z[i];

        tr_y_zz_xxz[i] = tr_y_0_xxz[i] * fe_0 + tr_y_z_xx[i] * fe_0 + tr_y_z_xxz[i] * pa_z[i];

        tr_y_zz_xyy[i] = tr_y_0_xyy[i] * fe_0 + tr_y_z_xyy[i] * pa_z[i];

        tr_y_zz_xyz[i] = tr_y_0_xyz[i] * fe_0 + tr_y_z_xy[i] * fe_0 + tr_y_z_xyz[i] * pa_z[i];

        tr_y_zz_xzz[i] = tr_y_0_xzz[i] * fe_0 + 2.0 * tr_y_z_xz[i] * fe_0 + tr_y_z_xzz[i] * pa_z[i];

        tr_y_zz_yyy[i] = tr_y_0_yyy[i] * fe_0 + tr_y_z_yyy[i] * pa_z[i];

        tr_y_zz_yyz[i] = tr_y_0_yyz[i] * fe_0 + tr_y_z_yy[i] * fe_0 + tr_y_z_yyz[i] * pa_z[i];

        tr_y_zz_yzz[i] = tr_y_0_yzz[i] * fe_0 + 2.0 * tr_y_z_yz[i] * fe_0 + tr_y_z_yzz[i] * pa_z[i];

        tr_y_zz_zzz[i] = tr_y_0_zzz[i] * fe_0 + 3.0 * tr_y_z_zz[i] * fe_0 + tr_y_z_zzz[i] * pa_z[i];
    }

    // Set up 120-130 components of targeted buffer : DF

    auto tr_z_xx_xxx = pbuffer.data(idx_dip_df + 120);

    auto tr_z_xx_xxy = pbuffer.data(idx_dip_df + 121);

    auto tr_z_xx_xxz = pbuffer.data(idx_dip_df + 122);

    auto tr_z_xx_xyy = pbuffer.data(idx_dip_df + 123);

    auto tr_z_xx_xyz = pbuffer.data(idx_dip_df + 124);

    auto tr_z_xx_xzz = pbuffer.data(idx_dip_df + 125);

    auto tr_z_xx_yyy = pbuffer.data(idx_dip_df + 126);

    auto tr_z_xx_yyz = pbuffer.data(idx_dip_df + 127);

    auto tr_z_xx_yzz = pbuffer.data(idx_dip_df + 128);

    auto tr_z_xx_zzz = pbuffer.data(idx_dip_df + 129);

#pragma omp simd aligned(pa_x,            \
                             tr_z_0_xxx,  \
                             tr_z_0_xxy,  \
                             tr_z_0_xxz,  \
                             tr_z_0_xyy,  \
                             tr_z_0_xyz,  \
                             tr_z_0_xzz,  \
                             tr_z_0_yyy,  \
                             tr_z_0_yyz,  \
                             tr_z_0_yzz,  \
                             tr_z_0_zzz,  \
                             tr_z_x_xx,   \
                             tr_z_x_xxx,  \
                             tr_z_x_xxy,  \
                             tr_z_x_xxz,  \
                             tr_z_x_xy,   \
                             tr_z_x_xyy,  \
                             tr_z_x_xyz,  \
                             tr_z_x_xz,   \
                             tr_z_x_xzz,  \
                             tr_z_x_yy,   \
                             tr_z_x_yyy,  \
                             tr_z_x_yyz,  \
                             tr_z_x_yz,   \
                             tr_z_x_yzz,  \
                             tr_z_x_zz,   \
                             tr_z_x_zzz,  \
                             tr_z_xx_xxx, \
                             tr_z_xx_xxy, \
                             tr_z_xx_xxz, \
                             tr_z_xx_xyy, \
                             tr_z_xx_xyz, \
                             tr_z_xx_xzz, \
                             tr_z_xx_yyy, \
                             tr_z_xx_yyz, \
                             tr_z_xx_yzz, \
                             tr_z_xx_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xx_xxx[i] = tr_z_0_xxx[i] * fe_0 + 3.0 * tr_z_x_xx[i] * fe_0 + tr_z_x_xxx[i] * pa_x[i];

        tr_z_xx_xxy[i] = tr_z_0_xxy[i] * fe_0 + 2.0 * tr_z_x_xy[i] * fe_0 + tr_z_x_xxy[i] * pa_x[i];

        tr_z_xx_xxz[i] = tr_z_0_xxz[i] * fe_0 + 2.0 * tr_z_x_xz[i] * fe_0 + tr_z_x_xxz[i] * pa_x[i];

        tr_z_xx_xyy[i] = tr_z_0_xyy[i] * fe_0 + tr_z_x_yy[i] * fe_0 + tr_z_x_xyy[i] * pa_x[i];

        tr_z_xx_xyz[i] = tr_z_0_xyz[i] * fe_0 + tr_z_x_yz[i] * fe_0 + tr_z_x_xyz[i] * pa_x[i];

        tr_z_xx_xzz[i] = tr_z_0_xzz[i] * fe_0 + tr_z_x_zz[i] * fe_0 + tr_z_x_xzz[i] * pa_x[i];

        tr_z_xx_yyy[i] = tr_z_0_yyy[i] * fe_0 + tr_z_x_yyy[i] * pa_x[i];

        tr_z_xx_yyz[i] = tr_z_0_yyz[i] * fe_0 + tr_z_x_yyz[i] * pa_x[i];

        tr_z_xx_yzz[i] = tr_z_0_yzz[i] * fe_0 + tr_z_x_yzz[i] * pa_x[i];

        tr_z_xx_zzz[i] = tr_z_0_zzz[i] * fe_0 + tr_z_x_zzz[i] * pa_x[i];
    }

    // Set up 130-140 components of targeted buffer : DF

    auto tr_z_xy_xxx = pbuffer.data(idx_dip_df + 130);

    auto tr_z_xy_xxy = pbuffer.data(idx_dip_df + 131);

    auto tr_z_xy_xxz = pbuffer.data(idx_dip_df + 132);

    auto tr_z_xy_xyy = pbuffer.data(idx_dip_df + 133);

    auto tr_z_xy_xyz = pbuffer.data(idx_dip_df + 134);

    auto tr_z_xy_xzz = pbuffer.data(idx_dip_df + 135);

    auto tr_z_xy_yyy = pbuffer.data(idx_dip_df + 136);

    auto tr_z_xy_yyz = pbuffer.data(idx_dip_df + 137);

    auto tr_z_xy_yzz = pbuffer.data(idx_dip_df + 138);

    auto tr_z_xy_zzz = pbuffer.data(idx_dip_df + 139);

#pragma omp simd aligned(pa_x,            \
                             pa_y,        \
                             tr_z_x_xxx,  \
                             tr_z_x_xxz,  \
                             tr_z_x_xzz,  \
                             tr_z_xy_xxx, \
                             tr_z_xy_xxy, \
                             tr_z_xy_xxz, \
                             tr_z_xy_xyy, \
                             tr_z_xy_xyz, \
                             tr_z_xy_xzz, \
                             tr_z_xy_yyy, \
                             tr_z_xy_yyz, \
                             tr_z_xy_yzz, \
                             tr_z_xy_zzz, \
                             tr_z_y_xxy,  \
                             tr_z_y_xy,   \
                             tr_z_y_xyy,  \
                             tr_z_y_xyz,  \
                             tr_z_y_yy,   \
                             tr_z_y_yyy,  \
                             tr_z_y_yyz,  \
                             tr_z_y_yz,   \
                             tr_z_y_yzz,  \
                             tr_z_y_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xy_xxx[i] = tr_z_x_xxx[i] * pa_y[i];

        tr_z_xy_xxy[i] = 2.0 * tr_z_y_xy[i] * fe_0 + tr_z_y_xxy[i] * pa_x[i];

        tr_z_xy_xxz[i] = tr_z_x_xxz[i] * pa_y[i];

        tr_z_xy_xyy[i] = tr_z_y_yy[i] * fe_0 + tr_z_y_xyy[i] * pa_x[i];

        tr_z_xy_xyz[i] = tr_z_y_yz[i] * fe_0 + tr_z_y_xyz[i] * pa_x[i];

        tr_z_xy_xzz[i] = tr_z_x_xzz[i] * pa_y[i];

        tr_z_xy_yyy[i] = tr_z_y_yyy[i] * pa_x[i];

        tr_z_xy_yyz[i] = tr_z_y_yyz[i] * pa_x[i];

        tr_z_xy_yzz[i] = tr_z_y_yzz[i] * pa_x[i];

        tr_z_xy_zzz[i] = tr_z_y_zzz[i] * pa_x[i];
    }

    // Set up 140-150 components of targeted buffer : DF

    auto tr_z_xz_xxx = pbuffer.data(idx_dip_df + 140);

    auto tr_z_xz_xxy = pbuffer.data(idx_dip_df + 141);

    auto tr_z_xz_xxz = pbuffer.data(idx_dip_df + 142);

    auto tr_z_xz_xyy = pbuffer.data(idx_dip_df + 143);

    auto tr_z_xz_xyz = pbuffer.data(idx_dip_df + 144);

    auto tr_z_xz_xzz = pbuffer.data(idx_dip_df + 145);

    auto tr_z_xz_yyy = pbuffer.data(idx_dip_df + 146);

    auto tr_z_xz_yyz = pbuffer.data(idx_dip_df + 147);

    auto tr_z_xz_yzz = pbuffer.data(idx_dip_df + 148);

    auto tr_z_xz_zzz = pbuffer.data(idx_dip_df + 149);

#pragma omp simd aligned(pa_x,            \
                             tr_z_xz_xxx, \
                             tr_z_xz_xxy, \
                             tr_z_xz_xxz, \
                             tr_z_xz_xyy, \
                             tr_z_xz_xyz, \
                             tr_z_xz_xzz, \
                             tr_z_xz_yyy, \
                             tr_z_xz_yyz, \
                             tr_z_xz_yzz, \
                             tr_z_xz_zzz, \
                             tr_z_z_xx,   \
                             tr_z_z_xxx,  \
                             tr_z_z_xxy,  \
                             tr_z_z_xxz,  \
                             tr_z_z_xy,   \
                             tr_z_z_xyy,  \
                             tr_z_z_xyz,  \
                             tr_z_z_xz,   \
                             tr_z_z_xzz,  \
                             tr_z_z_yy,   \
                             tr_z_z_yyy,  \
                             tr_z_z_yyz,  \
                             tr_z_z_yz,   \
                             tr_z_z_yzz,  \
                             tr_z_z_zz,   \
                             tr_z_z_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xz_xxx[i] = 3.0 * tr_z_z_xx[i] * fe_0 + tr_z_z_xxx[i] * pa_x[i];

        tr_z_xz_xxy[i] = 2.0 * tr_z_z_xy[i] * fe_0 + tr_z_z_xxy[i] * pa_x[i];

        tr_z_xz_xxz[i] = 2.0 * tr_z_z_xz[i] * fe_0 + tr_z_z_xxz[i] * pa_x[i];

        tr_z_xz_xyy[i] = tr_z_z_yy[i] * fe_0 + tr_z_z_xyy[i] * pa_x[i];

        tr_z_xz_xyz[i] = tr_z_z_yz[i] * fe_0 + tr_z_z_xyz[i] * pa_x[i];

        tr_z_xz_xzz[i] = tr_z_z_zz[i] * fe_0 + tr_z_z_xzz[i] * pa_x[i];

        tr_z_xz_yyy[i] = tr_z_z_yyy[i] * pa_x[i];

        tr_z_xz_yyz[i] = tr_z_z_yyz[i] * pa_x[i];

        tr_z_xz_yzz[i] = tr_z_z_yzz[i] * pa_x[i];

        tr_z_xz_zzz[i] = tr_z_z_zzz[i] * pa_x[i];
    }

    // Set up 150-160 components of targeted buffer : DF

    auto tr_z_yy_xxx = pbuffer.data(idx_dip_df + 150);

    auto tr_z_yy_xxy = pbuffer.data(idx_dip_df + 151);

    auto tr_z_yy_xxz = pbuffer.data(idx_dip_df + 152);

    auto tr_z_yy_xyy = pbuffer.data(idx_dip_df + 153);

    auto tr_z_yy_xyz = pbuffer.data(idx_dip_df + 154);

    auto tr_z_yy_xzz = pbuffer.data(idx_dip_df + 155);

    auto tr_z_yy_yyy = pbuffer.data(idx_dip_df + 156);

    auto tr_z_yy_yyz = pbuffer.data(idx_dip_df + 157);

    auto tr_z_yy_yzz = pbuffer.data(idx_dip_df + 158);

    auto tr_z_yy_zzz = pbuffer.data(idx_dip_df + 159);

#pragma omp simd aligned(pa_y,            \
                             tr_z_0_xxx,  \
                             tr_z_0_xxy,  \
                             tr_z_0_xxz,  \
                             tr_z_0_xyy,  \
                             tr_z_0_xyz,  \
                             tr_z_0_xzz,  \
                             tr_z_0_yyy,  \
                             tr_z_0_yyz,  \
                             tr_z_0_yzz,  \
                             tr_z_0_zzz,  \
                             tr_z_y_xx,   \
                             tr_z_y_xxx,  \
                             tr_z_y_xxy,  \
                             tr_z_y_xxz,  \
                             tr_z_y_xy,   \
                             tr_z_y_xyy,  \
                             tr_z_y_xyz,  \
                             tr_z_y_xz,   \
                             tr_z_y_xzz,  \
                             tr_z_y_yy,   \
                             tr_z_y_yyy,  \
                             tr_z_y_yyz,  \
                             tr_z_y_yz,   \
                             tr_z_y_yzz,  \
                             tr_z_y_zz,   \
                             tr_z_y_zzz,  \
                             tr_z_yy_xxx, \
                             tr_z_yy_xxy, \
                             tr_z_yy_xxz, \
                             tr_z_yy_xyy, \
                             tr_z_yy_xyz, \
                             tr_z_yy_xzz, \
                             tr_z_yy_yyy, \
                             tr_z_yy_yyz, \
                             tr_z_yy_yzz, \
                             tr_z_yy_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yy_xxx[i] = tr_z_0_xxx[i] * fe_0 + tr_z_y_xxx[i] * pa_y[i];

        tr_z_yy_xxy[i] = tr_z_0_xxy[i] * fe_0 + tr_z_y_xx[i] * fe_0 + tr_z_y_xxy[i] * pa_y[i];

        tr_z_yy_xxz[i] = tr_z_0_xxz[i] * fe_0 + tr_z_y_xxz[i] * pa_y[i];

        tr_z_yy_xyy[i] = tr_z_0_xyy[i] * fe_0 + 2.0 * tr_z_y_xy[i] * fe_0 + tr_z_y_xyy[i] * pa_y[i];

        tr_z_yy_xyz[i] = tr_z_0_xyz[i] * fe_0 + tr_z_y_xz[i] * fe_0 + tr_z_y_xyz[i] * pa_y[i];

        tr_z_yy_xzz[i] = tr_z_0_xzz[i] * fe_0 + tr_z_y_xzz[i] * pa_y[i];

        tr_z_yy_yyy[i] = tr_z_0_yyy[i] * fe_0 + 3.0 * tr_z_y_yy[i] * fe_0 + tr_z_y_yyy[i] * pa_y[i];

        tr_z_yy_yyz[i] = tr_z_0_yyz[i] * fe_0 + 2.0 * tr_z_y_yz[i] * fe_0 + tr_z_y_yyz[i] * pa_y[i];

        tr_z_yy_yzz[i] = tr_z_0_yzz[i] * fe_0 + tr_z_y_zz[i] * fe_0 + tr_z_y_yzz[i] * pa_y[i];

        tr_z_yy_zzz[i] = tr_z_0_zzz[i] * fe_0 + tr_z_y_zzz[i] * pa_y[i];
    }

    // Set up 160-170 components of targeted buffer : DF

    auto tr_z_yz_xxx = pbuffer.data(idx_dip_df + 160);

    auto tr_z_yz_xxy = pbuffer.data(idx_dip_df + 161);

    auto tr_z_yz_xxz = pbuffer.data(idx_dip_df + 162);

    auto tr_z_yz_xyy = pbuffer.data(idx_dip_df + 163);

    auto tr_z_yz_xyz = pbuffer.data(idx_dip_df + 164);

    auto tr_z_yz_xzz = pbuffer.data(idx_dip_df + 165);

    auto tr_z_yz_yyy = pbuffer.data(idx_dip_df + 166);

    auto tr_z_yz_yyz = pbuffer.data(idx_dip_df + 167);

    auto tr_z_yz_yzz = pbuffer.data(idx_dip_df + 168);

    auto tr_z_yz_zzz = pbuffer.data(idx_dip_df + 169);

#pragma omp simd aligned(pa_y,            \
                             tr_z_yz_xxx, \
                             tr_z_yz_xxy, \
                             tr_z_yz_xxz, \
                             tr_z_yz_xyy, \
                             tr_z_yz_xyz, \
                             tr_z_yz_xzz, \
                             tr_z_yz_yyy, \
                             tr_z_yz_yyz, \
                             tr_z_yz_yzz, \
                             tr_z_yz_zzz, \
                             tr_z_z_xx,   \
                             tr_z_z_xxx,  \
                             tr_z_z_xxy,  \
                             tr_z_z_xxz,  \
                             tr_z_z_xy,   \
                             tr_z_z_xyy,  \
                             tr_z_z_xyz,  \
                             tr_z_z_xz,   \
                             tr_z_z_xzz,  \
                             tr_z_z_yy,   \
                             tr_z_z_yyy,  \
                             tr_z_z_yyz,  \
                             tr_z_z_yz,   \
                             tr_z_z_yzz,  \
                             tr_z_z_zz,   \
                             tr_z_z_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yz_xxx[i] = tr_z_z_xxx[i] * pa_y[i];

        tr_z_yz_xxy[i] = tr_z_z_xx[i] * fe_0 + tr_z_z_xxy[i] * pa_y[i];

        tr_z_yz_xxz[i] = tr_z_z_xxz[i] * pa_y[i];

        tr_z_yz_xyy[i] = 2.0 * tr_z_z_xy[i] * fe_0 + tr_z_z_xyy[i] * pa_y[i];

        tr_z_yz_xyz[i] = tr_z_z_xz[i] * fe_0 + tr_z_z_xyz[i] * pa_y[i];

        tr_z_yz_xzz[i] = tr_z_z_xzz[i] * pa_y[i];

        tr_z_yz_yyy[i] = 3.0 * tr_z_z_yy[i] * fe_0 + tr_z_z_yyy[i] * pa_y[i];

        tr_z_yz_yyz[i] = 2.0 * tr_z_z_yz[i] * fe_0 + tr_z_z_yyz[i] * pa_y[i];

        tr_z_yz_yzz[i] = tr_z_z_zz[i] * fe_0 + tr_z_z_yzz[i] * pa_y[i];

        tr_z_yz_zzz[i] = tr_z_z_zzz[i] * pa_y[i];
    }

    // Set up 170-180 components of targeted buffer : DF

    auto tr_z_zz_xxx = pbuffer.data(idx_dip_df + 170);

    auto tr_z_zz_xxy = pbuffer.data(idx_dip_df + 171);

    auto tr_z_zz_xxz = pbuffer.data(idx_dip_df + 172);

    auto tr_z_zz_xyy = pbuffer.data(idx_dip_df + 173);

    auto tr_z_zz_xyz = pbuffer.data(idx_dip_df + 174);

    auto tr_z_zz_xzz = pbuffer.data(idx_dip_df + 175);

    auto tr_z_zz_yyy = pbuffer.data(idx_dip_df + 176);

    auto tr_z_zz_yyz = pbuffer.data(idx_dip_df + 177);

    auto tr_z_zz_yzz = pbuffer.data(idx_dip_df + 178);

    auto tr_z_zz_zzz = pbuffer.data(idx_dip_df + 179);

#pragma omp simd aligned(pa_z,            \
                             tr_z_0_xxx,  \
                             tr_z_0_xxy,  \
                             tr_z_0_xxz,  \
                             tr_z_0_xyy,  \
                             tr_z_0_xyz,  \
                             tr_z_0_xzz,  \
                             tr_z_0_yyy,  \
                             tr_z_0_yyz,  \
                             tr_z_0_yzz,  \
                             tr_z_0_zzz,  \
                             tr_z_z_xx,   \
                             tr_z_z_xxx,  \
                             tr_z_z_xxy,  \
                             tr_z_z_xxz,  \
                             tr_z_z_xy,   \
                             tr_z_z_xyy,  \
                             tr_z_z_xyz,  \
                             tr_z_z_xz,   \
                             tr_z_z_xzz,  \
                             tr_z_z_yy,   \
                             tr_z_z_yyy,  \
                             tr_z_z_yyz,  \
                             tr_z_z_yz,   \
                             tr_z_z_yzz,  \
                             tr_z_z_zz,   \
                             tr_z_z_zzz,  \
                             tr_z_zz_xxx, \
                             tr_z_zz_xxy, \
                             tr_z_zz_xxz, \
                             tr_z_zz_xyy, \
                             tr_z_zz_xyz, \
                             tr_z_zz_xzz, \
                             tr_z_zz_yyy, \
                             tr_z_zz_yyz, \
                             tr_z_zz_yzz, \
                             tr_z_zz_zzz, \
                             ts_z_xxx,    \
                             ts_z_xxy,    \
                             ts_z_xxz,    \
                             ts_z_xyy,    \
                             ts_z_xyz,    \
                             ts_z_xzz,    \
                             ts_z_yyy,    \
                             ts_z_yyz,    \
                             ts_z_yzz,    \
                             ts_z_zzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_zz_xxx[i] = tr_z_0_xxx[i] * fe_0 + ts_z_xxx[i] * fe_0 + tr_z_z_xxx[i] * pa_z[i];

        tr_z_zz_xxy[i] = tr_z_0_xxy[i] * fe_0 + ts_z_xxy[i] * fe_0 + tr_z_z_xxy[i] * pa_z[i];

        tr_z_zz_xxz[i] = tr_z_0_xxz[i] * fe_0 + tr_z_z_xx[i] * fe_0 + ts_z_xxz[i] * fe_0 + tr_z_z_xxz[i] * pa_z[i];

        tr_z_zz_xyy[i] = tr_z_0_xyy[i] * fe_0 + ts_z_xyy[i] * fe_0 + tr_z_z_xyy[i] * pa_z[i];

        tr_z_zz_xyz[i] = tr_z_0_xyz[i] * fe_0 + tr_z_z_xy[i] * fe_0 + ts_z_xyz[i] * fe_0 + tr_z_z_xyz[i] * pa_z[i];

        tr_z_zz_xzz[i] = tr_z_0_xzz[i] * fe_0 + 2.0 * tr_z_z_xz[i] * fe_0 + ts_z_xzz[i] * fe_0 + tr_z_z_xzz[i] * pa_z[i];

        tr_z_zz_yyy[i] = tr_z_0_yyy[i] * fe_0 + ts_z_yyy[i] * fe_0 + tr_z_z_yyy[i] * pa_z[i];

        tr_z_zz_yyz[i] = tr_z_0_yyz[i] * fe_0 + tr_z_z_yy[i] * fe_0 + ts_z_yyz[i] * fe_0 + tr_z_z_yyz[i] * pa_z[i];

        tr_z_zz_yzz[i] = tr_z_0_yzz[i] * fe_0 + 2.0 * tr_z_z_yz[i] * fe_0 + ts_z_yzz[i] * fe_0 + tr_z_z_yzz[i] * pa_z[i];

        tr_z_zz_zzz[i] = tr_z_0_zzz[i] * fe_0 + 3.0 * tr_z_z_zz[i] * fe_0 + ts_z_zzz[i] * fe_0 + tr_z_z_zzz[i] * pa_z[i];
    }
}

}  // namespace diprec
