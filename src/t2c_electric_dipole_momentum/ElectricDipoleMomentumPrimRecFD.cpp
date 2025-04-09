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

#include "ElectricDipoleMomentumPrimRecFD.hpp"

namespace diprec {  // diprec namespace

auto
comp_prim_electric_dipole_momentum_fd(CSimdArray<double>&       pbuffer,
                                      const size_t              idx_dip_fd,
                                      const size_t              idx_dip_pd,
                                      const size_t              idx_dip_dp,
                                      const size_t              idx_ovl_dd,
                                      const size_t              idx_dip_dd,
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

    // Set up components of auxiliary buffer : DP

    auto tr_x_xx_x = pbuffer.data(idx_dip_dp);

    auto tr_x_xx_y = pbuffer.data(idx_dip_dp + 1);

    auto tr_x_xx_z = pbuffer.data(idx_dip_dp + 2);

    auto tr_x_yy_x = pbuffer.data(idx_dip_dp + 9);

    auto tr_x_yy_y = pbuffer.data(idx_dip_dp + 10);

    auto tr_x_yy_z = pbuffer.data(idx_dip_dp + 11);

    auto tr_x_zz_x = pbuffer.data(idx_dip_dp + 15);

    auto tr_x_zz_y = pbuffer.data(idx_dip_dp + 16);

    auto tr_x_zz_z = pbuffer.data(idx_dip_dp + 17);

    auto tr_y_xx_x = pbuffer.data(idx_dip_dp + 18);

    auto tr_y_xx_y = pbuffer.data(idx_dip_dp + 19);

    auto tr_y_xx_z = pbuffer.data(idx_dip_dp + 20);

    auto tr_y_xy_y = pbuffer.data(idx_dip_dp + 22);

    auto tr_y_yy_x = pbuffer.data(idx_dip_dp + 27);

    auto tr_y_yy_y = pbuffer.data(idx_dip_dp + 28);

    auto tr_y_yy_z = pbuffer.data(idx_dip_dp + 29);

    auto tr_y_yz_z = pbuffer.data(idx_dip_dp + 32);

    auto tr_y_zz_x = pbuffer.data(idx_dip_dp + 33);

    auto tr_y_zz_y = pbuffer.data(idx_dip_dp + 34);

    auto tr_y_zz_z = pbuffer.data(idx_dip_dp + 35);

    auto tr_z_xx_x = pbuffer.data(idx_dip_dp + 36);

    auto tr_z_xx_y = pbuffer.data(idx_dip_dp + 37);

    auto tr_z_xx_z = pbuffer.data(idx_dip_dp + 38);

    auto tr_z_xz_z = pbuffer.data(idx_dip_dp + 44);

    auto tr_z_yy_x = pbuffer.data(idx_dip_dp + 45);

    auto tr_z_yy_y = pbuffer.data(idx_dip_dp + 46);

    auto tr_z_yy_z = pbuffer.data(idx_dip_dp + 47);

    auto tr_z_yz_y = pbuffer.data(idx_dip_dp + 49);

    auto tr_z_yz_z = pbuffer.data(idx_dip_dp + 50);

    auto tr_z_zz_x = pbuffer.data(idx_dip_dp + 51);

    auto tr_z_zz_y = pbuffer.data(idx_dip_dp + 52);

    auto tr_z_zz_z = pbuffer.data(idx_dip_dp + 53);

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

    // Set up components of auxiliary buffer : DD

    auto tr_x_xx_xx = pbuffer.data(idx_dip_dd);

    auto tr_x_xx_xy = pbuffer.data(idx_dip_dd + 1);

    auto tr_x_xx_xz = pbuffer.data(idx_dip_dd + 2);

    auto tr_x_xx_yy = pbuffer.data(idx_dip_dd + 3);

    auto tr_x_xx_yz = pbuffer.data(idx_dip_dd + 4);

    auto tr_x_xx_zz = pbuffer.data(idx_dip_dd + 5);

    auto tr_x_xy_xx = pbuffer.data(idx_dip_dd + 6);

    auto tr_x_xy_xy = pbuffer.data(idx_dip_dd + 7);

    auto tr_x_xy_xz = pbuffer.data(idx_dip_dd + 8);

    auto tr_x_xy_yy = pbuffer.data(idx_dip_dd + 9);

    auto tr_x_xz_xx = pbuffer.data(idx_dip_dd + 12);

    auto tr_x_xz_xy = pbuffer.data(idx_dip_dd + 13);

    auto tr_x_xz_xz = pbuffer.data(idx_dip_dd + 14);

    auto tr_x_xz_zz = pbuffer.data(idx_dip_dd + 17);

    auto tr_x_yy_xx = pbuffer.data(idx_dip_dd + 18);

    auto tr_x_yy_xy = pbuffer.data(idx_dip_dd + 19);

    auto tr_x_yy_xz = pbuffer.data(idx_dip_dd + 20);

    auto tr_x_yy_yy = pbuffer.data(idx_dip_dd + 21);

    auto tr_x_yy_yz = pbuffer.data(idx_dip_dd + 22);

    auto tr_x_yy_zz = pbuffer.data(idx_dip_dd + 23);

    auto tr_x_yz_xz = pbuffer.data(idx_dip_dd + 26);

    auto tr_x_yz_yz = pbuffer.data(idx_dip_dd + 28);

    auto tr_x_yz_zz = pbuffer.data(idx_dip_dd + 29);

    auto tr_x_zz_xx = pbuffer.data(idx_dip_dd + 30);

    auto tr_x_zz_xy = pbuffer.data(idx_dip_dd + 31);

    auto tr_x_zz_xz = pbuffer.data(idx_dip_dd + 32);

    auto tr_x_zz_yy = pbuffer.data(idx_dip_dd + 33);

    auto tr_x_zz_yz = pbuffer.data(idx_dip_dd + 34);

    auto tr_x_zz_zz = pbuffer.data(idx_dip_dd + 35);

    auto tr_y_xx_xx = pbuffer.data(idx_dip_dd + 36);

    auto tr_y_xx_xy = pbuffer.data(idx_dip_dd + 37);

    auto tr_y_xx_xz = pbuffer.data(idx_dip_dd + 38);

    auto tr_y_xx_yy = pbuffer.data(idx_dip_dd + 39);

    auto tr_y_xx_yz = pbuffer.data(idx_dip_dd + 40);

    auto tr_y_xx_zz = pbuffer.data(idx_dip_dd + 41);

    auto tr_y_xy_xx = pbuffer.data(idx_dip_dd + 42);

    auto tr_y_xy_xy = pbuffer.data(idx_dip_dd + 43);

    auto tr_y_xy_yy = pbuffer.data(idx_dip_dd + 45);

    auto tr_y_xy_yz = pbuffer.data(idx_dip_dd + 46);

    auto tr_y_xy_zz = pbuffer.data(idx_dip_dd + 47);

    auto tr_y_xz_yz = pbuffer.data(idx_dip_dd + 52);

    auto tr_y_xz_zz = pbuffer.data(idx_dip_dd + 53);

    auto tr_y_yy_xx = pbuffer.data(idx_dip_dd + 54);

    auto tr_y_yy_xy = pbuffer.data(idx_dip_dd + 55);

    auto tr_y_yy_xz = pbuffer.data(idx_dip_dd + 56);

    auto tr_y_yy_yy = pbuffer.data(idx_dip_dd + 57);

    auto tr_y_yy_yz = pbuffer.data(idx_dip_dd + 58);

    auto tr_y_yy_zz = pbuffer.data(idx_dip_dd + 59);

    auto tr_y_yz_xy = pbuffer.data(idx_dip_dd + 61);

    auto tr_y_yz_xz = pbuffer.data(idx_dip_dd + 62);

    auto tr_y_yz_yy = pbuffer.data(idx_dip_dd + 63);

    auto tr_y_yz_yz = pbuffer.data(idx_dip_dd + 64);

    auto tr_y_yz_zz = pbuffer.data(idx_dip_dd + 65);

    auto tr_y_zz_xx = pbuffer.data(idx_dip_dd + 66);

    auto tr_y_zz_xy = pbuffer.data(idx_dip_dd + 67);

    auto tr_y_zz_xz = pbuffer.data(idx_dip_dd + 68);

    auto tr_y_zz_yy = pbuffer.data(idx_dip_dd + 69);

    auto tr_y_zz_yz = pbuffer.data(idx_dip_dd + 70);

    auto tr_y_zz_zz = pbuffer.data(idx_dip_dd + 71);

    auto tr_z_xx_xx = pbuffer.data(idx_dip_dd + 72);

    auto tr_z_xx_xy = pbuffer.data(idx_dip_dd + 73);

    auto tr_z_xx_xz = pbuffer.data(idx_dip_dd + 74);

    auto tr_z_xx_yy = pbuffer.data(idx_dip_dd + 75);

    auto tr_z_xx_yz = pbuffer.data(idx_dip_dd + 76);

    auto tr_z_xx_zz = pbuffer.data(idx_dip_dd + 77);

    auto tr_z_xy_yy = pbuffer.data(idx_dip_dd + 81);

    auto tr_z_xy_yz = pbuffer.data(idx_dip_dd + 82);

    auto tr_z_xz_xx = pbuffer.data(idx_dip_dd + 84);

    auto tr_z_xz_xz = pbuffer.data(idx_dip_dd + 86);

    auto tr_z_xz_yy = pbuffer.data(idx_dip_dd + 87);

    auto tr_z_xz_yz = pbuffer.data(idx_dip_dd + 88);

    auto tr_z_xz_zz = pbuffer.data(idx_dip_dd + 89);

    auto tr_z_yy_xx = pbuffer.data(idx_dip_dd + 90);

    auto tr_z_yy_xy = pbuffer.data(idx_dip_dd + 91);

    auto tr_z_yy_xz = pbuffer.data(idx_dip_dd + 92);

    auto tr_z_yy_yy = pbuffer.data(idx_dip_dd + 93);

    auto tr_z_yy_yz = pbuffer.data(idx_dip_dd + 94);

    auto tr_z_yy_zz = pbuffer.data(idx_dip_dd + 95);

    auto tr_z_yz_xx = pbuffer.data(idx_dip_dd + 96);

    auto tr_z_yz_xy = pbuffer.data(idx_dip_dd + 97);

    auto tr_z_yz_xz = pbuffer.data(idx_dip_dd + 98);

    auto tr_z_yz_yy = pbuffer.data(idx_dip_dd + 99);

    auto tr_z_yz_yz = pbuffer.data(idx_dip_dd + 100);

    auto tr_z_yz_zz = pbuffer.data(idx_dip_dd + 101);

    auto tr_z_zz_xx = pbuffer.data(idx_dip_dd + 102);

    auto tr_z_zz_xy = pbuffer.data(idx_dip_dd + 103);

    auto tr_z_zz_xz = pbuffer.data(idx_dip_dd + 104);

    auto tr_z_zz_yy = pbuffer.data(idx_dip_dd + 105);

    auto tr_z_zz_yz = pbuffer.data(idx_dip_dd + 106);

    auto tr_z_zz_zz = pbuffer.data(idx_dip_dd + 107);

    // Set up 0-6 components of targeted buffer : FD

    auto tr_x_xxx_xx = pbuffer.data(idx_dip_fd);

    auto tr_x_xxx_xy = pbuffer.data(idx_dip_fd + 1);

    auto tr_x_xxx_xz = pbuffer.data(idx_dip_fd + 2);

    auto tr_x_xxx_yy = pbuffer.data(idx_dip_fd + 3);

    auto tr_x_xxx_yz = pbuffer.data(idx_dip_fd + 4);

    auto tr_x_xxx_zz = pbuffer.data(idx_dip_fd + 5);

#pragma omp simd aligned(pa_x,            \
                             tr_x_x_xx,   \
                             tr_x_x_xy,   \
                             tr_x_x_xz,   \
                             tr_x_x_yy,   \
                             tr_x_x_yz,   \
                             tr_x_x_zz,   \
                             tr_x_xx_x,   \
                             tr_x_xx_xx,  \
                             tr_x_xx_xy,  \
                             tr_x_xx_xz,  \
                             tr_x_xx_y,   \
                             tr_x_xx_yy,  \
                             tr_x_xx_yz,  \
                             tr_x_xx_z,   \
                             tr_x_xx_zz,  \
                             tr_x_xxx_xx, \
                             tr_x_xxx_xy, \
                             tr_x_xxx_xz, \
                             tr_x_xxx_yy, \
                             tr_x_xxx_yz, \
                             tr_x_xxx_zz, \
                             ts_xx_xx,    \
                             ts_xx_xy,    \
                             ts_xx_xz,    \
                             ts_xx_yy,    \
                             ts_xx_yz,    \
                             ts_xx_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxx_xx[i] = 2.0 * tr_x_x_xx[i] * fe_0 + 2.0 * tr_x_xx_x[i] * fe_0 + ts_xx_xx[i] * fe_0 + tr_x_xx_xx[i] * pa_x[i];

        tr_x_xxx_xy[i] = 2.0 * tr_x_x_xy[i] * fe_0 + tr_x_xx_y[i] * fe_0 + ts_xx_xy[i] * fe_0 + tr_x_xx_xy[i] * pa_x[i];

        tr_x_xxx_xz[i] = 2.0 * tr_x_x_xz[i] * fe_0 + tr_x_xx_z[i] * fe_0 + ts_xx_xz[i] * fe_0 + tr_x_xx_xz[i] * pa_x[i];

        tr_x_xxx_yy[i] = 2.0 * tr_x_x_yy[i] * fe_0 + ts_xx_yy[i] * fe_0 + tr_x_xx_yy[i] * pa_x[i];

        tr_x_xxx_yz[i] = 2.0 * tr_x_x_yz[i] * fe_0 + ts_xx_yz[i] * fe_0 + tr_x_xx_yz[i] * pa_x[i];

        tr_x_xxx_zz[i] = 2.0 * tr_x_x_zz[i] * fe_0 + ts_xx_zz[i] * fe_0 + tr_x_xx_zz[i] * pa_x[i];
    }

    // Set up 6-12 components of targeted buffer : FD

    auto tr_x_xxy_xx = pbuffer.data(idx_dip_fd + 6);

    auto tr_x_xxy_xy = pbuffer.data(idx_dip_fd + 7);

    auto tr_x_xxy_xz = pbuffer.data(idx_dip_fd + 8);

    auto tr_x_xxy_yy = pbuffer.data(idx_dip_fd + 9);

    auto tr_x_xxy_yz = pbuffer.data(idx_dip_fd + 10);

    auto tr_x_xxy_zz = pbuffer.data(idx_dip_fd + 11);

#pragma omp simd aligned(pa_y,            \
                             tr_x_xx_x,   \
                             tr_x_xx_xx,  \
                             tr_x_xx_xy,  \
                             tr_x_xx_xz,  \
                             tr_x_xx_y,   \
                             tr_x_xx_yy,  \
                             tr_x_xx_yz,  \
                             tr_x_xx_z,   \
                             tr_x_xx_zz,  \
                             tr_x_xxy_xx, \
                             tr_x_xxy_xy, \
                             tr_x_xxy_xz, \
                             tr_x_xxy_yy, \
                             tr_x_xxy_yz, \
                             tr_x_xxy_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxy_xx[i] = tr_x_xx_xx[i] * pa_y[i];

        tr_x_xxy_xy[i] = tr_x_xx_x[i] * fe_0 + tr_x_xx_xy[i] * pa_y[i];

        tr_x_xxy_xz[i] = tr_x_xx_xz[i] * pa_y[i];

        tr_x_xxy_yy[i] = 2.0 * tr_x_xx_y[i] * fe_0 + tr_x_xx_yy[i] * pa_y[i];

        tr_x_xxy_yz[i] = tr_x_xx_z[i] * fe_0 + tr_x_xx_yz[i] * pa_y[i];

        tr_x_xxy_zz[i] = tr_x_xx_zz[i] * pa_y[i];
    }

    // Set up 12-18 components of targeted buffer : FD

    auto tr_x_xxz_xx = pbuffer.data(idx_dip_fd + 12);

    auto tr_x_xxz_xy = pbuffer.data(idx_dip_fd + 13);

    auto tr_x_xxz_xz = pbuffer.data(idx_dip_fd + 14);

    auto tr_x_xxz_yy = pbuffer.data(idx_dip_fd + 15);

    auto tr_x_xxz_yz = pbuffer.data(idx_dip_fd + 16);

    auto tr_x_xxz_zz = pbuffer.data(idx_dip_fd + 17);

#pragma omp simd aligned(pa_z,            \
                             tr_x_xx_x,   \
                             tr_x_xx_xx,  \
                             tr_x_xx_xy,  \
                             tr_x_xx_xz,  \
                             tr_x_xx_y,   \
                             tr_x_xx_yy,  \
                             tr_x_xx_yz,  \
                             tr_x_xx_z,   \
                             tr_x_xx_zz,  \
                             tr_x_xxz_xx, \
                             tr_x_xxz_xy, \
                             tr_x_xxz_xz, \
                             tr_x_xxz_yy, \
                             tr_x_xxz_yz, \
                             tr_x_xxz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxz_xx[i] = tr_x_xx_xx[i] * pa_z[i];

        tr_x_xxz_xy[i] = tr_x_xx_xy[i] * pa_z[i];

        tr_x_xxz_xz[i] = tr_x_xx_x[i] * fe_0 + tr_x_xx_xz[i] * pa_z[i];

        tr_x_xxz_yy[i] = tr_x_xx_yy[i] * pa_z[i];

        tr_x_xxz_yz[i] = tr_x_xx_y[i] * fe_0 + tr_x_xx_yz[i] * pa_z[i];

        tr_x_xxz_zz[i] = 2.0 * tr_x_xx_z[i] * fe_0 + tr_x_xx_zz[i] * pa_z[i];
    }

    // Set up 18-24 components of targeted buffer : FD

    auto tr_x_xyy_xx = pbuffer.data(idx_dip_fd + 18);

    auto tr_x_xyy_xy = pbuffer.data(idx_dip_fd + 19);

    auto tr_x_xyy_xz = pbuffer.data(idx_dip_fd + 20);

    auto tr_x_xyy_yy = pbuffer.data(idx_dip_fd + 21);

    auto tr_x_xyy_yz = pbuffer.data(idx_dip_fd + 22);

    auto tr_x_xyy_zz = pbuffer.data(idx_dip_fd + 23);

#pragma omp simd aligned(pa_x,            \
                             pa_y,        \
                             tr_x_x_xx,   \
                             tr_x_x_xz,   \
                             tr_x_xy_xx,  \
                             tr_x_xy_xz,  \
                             tr_x_xyy_xx, \
                             tr_x_xyy_xy, \
                             tr_x_xyy_xz, \
                             tr_x_xyy_yy, \
                             tr_x_xyy_yz, \
                             tr_x_xyy_zz, \
                             tr_x_yy_xy,  \
                             tr_x_yy_y,   \
                             tr_x_yy_yy,  \
                             tr_x_yy_yz,  \
                             tr_x_yy_zz,  \
                             ts_yy_xy,    \
                             ts_yy_yy,    \
                             ts_yy_yz,    \
                             ts_yy_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyy_xx[i] = tr_x_x_xx[i] * fe_0 + tr_x_xy_xx[i] * pa_y[i];

        tr_x_xyy_xy[i] = tr_x_yy_y[i] * fe_0 + ts_yy_xy[i] * fe_0 + tr_x_yy_xy[i] * pa_x[i];

        tr_x_xyy_xz[i] = tr_x_x_xz[i] * fe_0 + tr_x_xy_xz[i] * pa_y[i];

        tr_x_xyy_yy[i] = ts_yy_yy[i] * fe_0 + tr_x_yy_yy[i] * pa_x[i];

        tr_x_xyy_yz[i] = ts_yy_yz[i] * fe_0 + tr_x_yy_yz[i] * pa_x[i];

        tr_x_xyy_zz[i] = ts_yy_zz[i] * fe_0 + tr_x_yy_zz[i] * pa_x[i];
    }

    // Set up 24-30 components of targeted buffer : FD

    auto tr_x_xyz_xx = pbuffer.data(idx_dip_fd + 24);

    auto tr_x_xyz_xy = pbuffer.data(idx_dip_fd + 25);

    auto tr_x_xyz_xz = pbuffer.data(idx_dip_fd + 26);

    auto tr_x_xyz_yy = pbuffer.data(idx_dip_fd + 27);

    auto tr_x_xyz_yz = pbuffer.data(idx_dip_fd + 28);

    auto tr_x_xyz_zz = pbuffer.data(idx_dip_fd + 29);

#pragma omp simd aligned(pa_x,            \
                             pa_y,        \
                             pa_z,        \
                             tr_x_xy_xy,  \
                             tr_x_xy_yy,  \
                             tr_x_xyz_xx, \
                             tr_x_xyz_xy, \
                             tr_x_xyz_xz, \
                             tr_x_xyz_yy, \
                             tr_x_xyz_yz, \
                             tr_x_xyz_zz, \
                             tr_x_xz_xx,  \
                             tr_x_xz_xz,  \
                             tr_x_xz_zz,  \
                             tr_x_yz_yz,  \
                             ts_yz_yz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyz_xx[i] = tr_x_xz_xx[i] * pa_y[i];

        tr_x_xyz_xy[i] = tr_x_xy_xy[i] * pa_z[i];

        tr_x_xyz_xz[i] = tr_x_xz_xz[i] * pa_y[i];

        tr_x_xyz_yy[i] = tr_x_xy_yy[i] * pa_z[i];

        tr_x_xyz_yz[i] = ts_yz_yz[i] * fe_0 + tr_x_yz_yz[i] * pa_x[i];

        tr_x_xyz_zz[i] = tr_x_xz_zz[i] * pa_y[i];
    }

    // Set up 30-36 components of targeted buffer : FD

    auto tr_x_xzz_xx = pbuffer.data(idx_dip_fd + 30);

    auto tr_x_xzz_xy = pbuffer.data(idx_dip_fd + 31);

    auto tr_x_xzz_xz = pbuffer.data(idx_dip_fd + 32);

    auto tr_x_xzz_yy = pbuffer.data(idx_dip_fd + 33);

    auto tr_x_xzz_yz = pbuffer.data(idx_dip_fd + 34);

    auto tr_x_xzz_zz = pbuffer.data(idx_dip_fd + 35);

#pragma omp simd aligned(pa_x,            \
                             pa_z,        \
                             tr_x_x_xx,   \
                             tr_x_x_xy,   \
                             tr_x_xz_xx,  \
                             tr_x_xz_xy,  \
                             tr_x_xzz_xx, \
                             tr_x_xzz_xy, \
                             tr_x_xzz_xz, \
                             tr_x_xzz_yy, \
                             tr_x_xzz_yz, \
                             tr_x_xzz_zz, \
                             tr_x_zz_xz,  \
                             tr_x_zz_yy,  \
                             tr_x_zz_yz,  \
                             tr_x_zz_z,   \
                             tr_x_zz_zz,  \
                             ts_zz_xz,    \
                             ts_zz_yy,    \
                             ts_zz_yz,    \
                             ts_zz_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xzz_xx[i] = tr_x_x_xx[i] * fe_0 + tr_x_xz_xx[i] * pa_z[i];

        tr_x_xzz_xy[i] = tr_x_x_xy[i] * fe_0 + tr_x_xz_xy[i] * pa_z[i];

        tr_x_xzz_xz[i] = tr_x_zz_z[i] * fe_0 + ts_zz_xz[i] * fe_0 + tr_x_zz_xz[i] * pa_x[i];

        tr_x_xzz_yy[i] = ts_zz_yy[i] * fe_0 + tr_x_zz_yy[i] * pa_x[i];

        tr_x_xzz_yz[i] = ts_zz_yz[i] * fe_0 + tr_x_zz_yz[i] * pa_x[i];

        tr_x_xzz_zz[i] = ts_zz_zz[i] * fe_0 + tr_x_zz_zz[i] * pa_x[i];
    }

    // Set up 36-42 components of targeted buffer : FD

    auto tr_x_yyy_xx = pbuffer.data(idx_dip_fd + 36);

    auto tr_x_yyy_xy = pbuffer.data(idx_dip_fd + 37);

    auto tr_x_yyy_xz = pbuffer.data(idx_dip_fd + 38);

    auto tr_x_yyy_yy = pbuffer.data(idx_dip_fd + 39);

    auto tr_x_yyy_yz = pbuffer.data(idx_dip_fd + 40);

    auto tr_x_yyy_zz = pbuffer.data(idx_dip_fd + 41);

#pragma omp simd aligned(pa_y,            \
                             tr_x_y_xx,   \
                             tr_x_y_xy,   \
                             tr_x_y_xz,   \
                             tr_x_y_yy,   \
                             tr_x_y_yz,   \
                             tr_x_y_zz,   \
                             tr_x_yy_x,   \
                             tr_x_yy_xx,  \
                             tr_x_yy_xy,  \
                             tr_x_yy_xz,  \
                             tr_x_yy_y,   \
                             tr_x_yy_yy,  \
                             tr_x_yy_yz,  \
                             tr_x_yy_z,   \
                             tr_x_yy_zz,  \
                             tr_x_yyy_xx, \
                             tr_x_yyy_xy, \
                             tr_x_yyy_xz, \
                             tr_x_yyy_yy, \
                             tr_x_yyy_yz, \
                             tr_x_yyy_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyy_xx[i] = 2.0 * tr_x_y_xx[i] * fe_0 + tr_x_yy_xx[i] * pa_y[i];

        tr_x_yyy_xy[i] = 2.0 * tr_x_y_xy[i] * fe_0 + tr_x_yy_x[i] * fe_0 + tr_x_yy_xy[i] * pa_y[i];

        tr_x_yyy_xz[i] = 2.0 * tr_x_y_xz[i] * fe_0 + tr_x_yy_xz[i] * pa_y[i];

        tr_x_yyy_yy[i] = 2.0 * tr_x_y_yy[i] * fe_0 + 2.0 * tr_x_yy_y[i] * fe_0 + tr_x_yy_yy[i] * pa_y[i];

        tr_x_yyy_yz[i] = 2.0 * tr_x_y_yz[i] * fe_0 + tr_x_yy_z[i] * fe_0 + tr_x_yy_yz[i] * pa_y[i];

        tr_x_yyy_zz[i] = 2.0 * tr_x_y_zz[i] * fe_0 + tr_x_yy_zz[i] * pa_y[i];
    }

    // Set up 42-48 components of targeted buffer : FD

    auto tr_x_yyz_xx = pbuffer.data(idx_dip_fd + 42);

    auto tr_x_yyz_xy = pbuffer.data(idx_dip_fd + 43);

    auto tr_x_yyz_xz = pbuffer.data(idx_dip_fd + 44);

    auto tr_x_yyz_yy = pbuffer.data(idx_dip_fd + 45);

    auto tr_x_yyz_yz = pbuffer.data(idx_dip_fd + 46);

    auto tr_x_yyz_zz = pbuffer.data(idx_dip_fd + 47);

#pragma omp simd aligned(pa_y,            \
                             pa_z,        \
                             tr_x_yy_xx,  \
                             tr_x_yy_xy,  \
                             tr_x_yy_y,   \
                             tr_x_yy_yy,  \
                             tr_x_yy_yz,  \
                             tr_x_yyz_xx, \
                             tr_x_yyz_xy, \
                             tr_x_yyz_xz, \
                             tr_x_yyz_yy, \
                             tr_x_yyz_yz, \
                             tr_x_yyz_zz, \
                             tr_x_yz_xz,  \
                             tr_x_yz_zz,  \
                             tr_x_z_xz,   \
                             tr_x_z_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyz_xx[i] = tr_x_yy_xx[i] * pa_z[i];

        tr_x_yyz_xy[i] = tr_x_yy_xy[i] * pa_z[i];

        tr_x_yyz_xz[i] = tr_x_z_xz[i] * fe_0 + tr_x_yz_xz[i] * pa_y[i];

        tr_x_yyz_yy[i] = tr_x_yy_yy[i] * pa_z[i];

        tr_x_yyz_yz[i] = tr_x_yy_y[i] * fe_0 + tr_x_yy_yz[i] * pa_z[i];

        tr_x_yyz_zz[i] = tr_x_z_zz[i] * fe_0 + tr_x_yz_zz[i] * pa_y[i];
    }

    // Set up 48-54 components of targeted buffer : FD

    auto tr_x_yzz_xx = pbuffer.data(idx_dip_fd + 48);

    auto tr_x_yzz_xy = pbuffer.data(idx_dip_fd + 49);

    auto tr_x_yzz_xz = pbuffer.data(idx_dip_fd + 50);

    auto tr_x_yzz_yy = pbuffer.data(idx_dip_fd + 51);

    auto tr_x_yzz_yz = pbuffer.data(idx_dip_fd + 52);

    auto tr_x_yzz_zz = pbuffer.data(idx_dip_fd + 53);

#pragma omp simd aligned(pa_y,            \
                             tr_x_yzz_xx, \
                             tr_x_yzz_xy, \
                             tr_x_yzz_xz, \
                             tr_x_yzz_yy, \
                             tr_x_yzz_yz, \
                             tr_x_yzz_zz, \
                             tr_x_zz_x,   \
                             tr_x_zz_xx,  \
                             tr_x_zz_xy,  \
                             tr_x_zz_xz,  \
                             tr_x_zz_y,   \
                             tr_x_zz_yy,  \
                             tr_x_zz_yz,  \
                             tr_x_zz_z,   \
                             tr_x_zz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yzz_xx[i] = tr_x_zz_xx[i] * pa_y[i];

        tr_x_yzz_xy[i] = tr_x_zz_x[i] * fe_0 + tr_x_zz_xy[i] * pa_y[i];

        tr_x_yzz_xz[i] = tr_x_zz_xz[i] * pa_y[i];

        tr_x_yzz_yy[i] = 2.0 * tr_x_zz_y[i] * fe_0 + tr_x_zz_yy[i] * pa_y[i];

        tr_x_yzz_yz[i] = tr_x_zz_z[i] * fe_0 + tr_x_zz_yz[i] * pa_y[i];

        tr_x_yzz_zz[i] = tr_x_zz_zz[i] * pa_y[i];
    }

    // Set up 54-60 components of targeted buffer : FD

    auto tr_x_zzz_xx = pbuffer.data(idx_dip_fd + 54);

    auto tr_x_zzz_xy = pbuffer.data(idx_dip_fd + 55);

    auto tr_x_zzz_xz = pbuffer.data(idx_dip_fd + 56);

    auto tr_x_zzz_yy = pbuffer.data(idx_dip_fd + 57);

    auto tr_x_zzz_yz = pbuffer.data(idx_dip_fd + 58);

    auto tr_x_zzz_zz = pbuffer.data(idx_dip_fd + 59);

#pragma omp simd aligned(pa_z,            \
                             tr_x_z_xx,   \
                             tr_x_z_xy,   \
                             tr_x_z_xz,   \
                             tr_x_z_yy,   \
                             tr_x_z_yz,   \
                             tr_x_z_zz,   \
                             tr_x_zz_x,   \
                             tr_x_zz_xx,  \
                             tr_x_zz_xy,  \
                             tr_x_zz_xz,  \
                             tr_x_zz_y,   \
                             tr_x_zz_yy,  \
                             tr_x_zz_yz,  \
                             tr_x_zz_z,   \
                             tr_x_zz_zz,  \
                             tr_x_zzz_xx, \
                             tr_x_zzz_xy, \
                             tr_x_zzz_xz, \
                             tr_x_zzz_yy, \
                             tr_x_zzz_yz, \
                             tr_x_zzz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_zzz_xx[i] = 2.0 * tr_x_z_xx[i] * fe_0 + tr_x_zz_xx[i] * pa_z[i];

        tr_x_zzz_xy[i] = 2.0 * tr_x_z_xy[i] * fe_0 + tr_x_zz_xy[i] * pa_z[i];

        tr_x_zzz_xz[i] = 2.0 * tr_x_z_xz[i] * fe_0 + tr_x_zz_x[i] * fe_0 + tr_x_zz_xz[i] * pa_z[i];

        tr_x_zzz_yy[i] = 2.0 * tr_x_z_yy[i] * fe_0 + tr_x_zz_yy[i] * pa_z[i];

        tr_x_zzz_yz[i] = 2.0 * tr_x_z_yz[i] * fe_0 + tr_x_zz_y[i] * fe_0 + tr_x_zz_yz[i] * pa_z[i];

        tr_x_zzz_zz[i] = 2.0 * tr_x_z_zz[i] * fe_0 + 2.0 * tr_x_zz_z[i] * fe_0 + tr_x_zz_zz[i] * pa_z[i];
    }

    // Set up 60-66 components of targeted buffer : FD

    auto tr_y_xxx_xx = pbuffer.data(idx_dip_fd + 60);

    auto tr_y_xxx_xy = pbuffer.data(idx_dip_fd + 61);

    auto tr_y_xxx_xz = pbuffer.data(idx_dip_fd + 62);

    auto tr_y_xxx_yy = pbuffer.data(idx_dip_fd + 63);

    auto tr_y_xxx_yz = pbuffer.data(idx_dip_fd + 64);

    auto tr_y_xxx_zz = pbuffer.data(idx_dip_fd + 65);

#pragma omp simd aligned(pa_x,            \
                             tr_y_x_xx,   \
                             tr_y_x_xy,   \
                             tr_y_x_xz,   \
                             tr_y_x_yy,   \
                             tr_y_x_yz,   \
                             tr_y_x_zz,   \
                             tr_y_xx_x,   \
                             tr_y_xx_xx,  \
                             tr_y_xx_xy,  \
                             tr_y_xx_xz,  \
                             tr_y_xx_y,   \
                             tr_y_xx_yy,  \
                             tr_y_xx_yz,  \
                             tr_y_xx_z,   \
                             tr_y_xx_zz,  \
                             tr_y_xxx_xx, \
                             tr_y_xxx_xy, \
                             tr_y_xxx_xz, \
                             tr_y_xxx_yy, \
                             tr_y_xxx_yz, \
                             tr_y_xxx_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxx_xx[i] = 2.0 * tr_y_x_xx[i] * fe_0 + 2.0 * tr_y_xx_x[i] * fe_0 + tr_y_xx_xx[i] * pa_x[i];

        tr_y_xxx_xy[i] = 2.0 * tr_y_x_xy[i] * fe_0 + tr_y_xx_y[i] * fe_0 + tr_y_xx_xy[i] * pa_x[i];

        tr_y_xxx_xz[i] = 2.0 * tr_y_x_xz[i] * fe_0 + tr_y_xx_z[i] * fe_0 + tr_y_xx_xz[i] * pa_x[i];

        tr_y_xxx_yy[i] = 2.0 * tr_y_x_yy[i] * fe_0 + tr_y_xx_yy[i] * pa_x[i];

        tr_y_xxx_yz[i] = 2.0 * tr_y_x_yz[i] * fe_0 + tr_y_xx_yz[i] * pa_x[i];

        tr_y_xxx_zz[i] = 2.0 * tr_y_x_zz[i] * fe_0 + tr_y_xx_zz[i] * pa_x[i];
    }

    // Set up 66-72 components of targeted buffer : FD

    auto tr_y_xxy_xx = pbuffer.data(idx_dip_fd + 66);

    auto tr_y_xxy_xy = pbuffer.data(idx_dip_fd + 67);

    auto tr_y_xxy_xz = pbuffer.data(idx_dip_fd + 68);

    auto tr_y_xxy_yy = pbuffer.data(idx_dip_fd + 69);

    auto tr_y_xxy_yz = pbuffer.data(idx_dip_fd + 70);

    auto tr_y_xxy_zz = pbuffer.data(idx_dip_fd + 71);

#pragma omp simd aligned(pa_x,            \
                             pa_y,        \
                             tr_y_xx_xx,  \
                             tr_y_xx_xz,  \
                             tr_y_xxy_xx, \
                             tr_y_xxy_xy, \
                             tr_y_xxy_xz, \
                             tr_y_xxy_yy, \
                             tr_y_xxy_yz, \
                             tr_y_xxy_zz, \
                             tr_y_xy_xy,  \
                             tr_y_xy_y,   \
                             tr_y_xy_yy,  \
                             tr_y_xy_yz,  \
                             tr_y_xy_zz,  \
                             tr_y_y_xy,   \
                             tr_y_y_yy,   \
                             tr_y_y_yz,   \
                             tr_y_y_zz,   \
                             ts_xx_xx,    \
                             ts_xx_xz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxy_xx[i] = ts_xx_xx[i] * fe_0 + tr_y_xx_xx[i] * pa_y[i];

        tr_y_xxy_xy[i] = tr_y_y_xy[i] * fe_0 + tr_y_xy_y[i] * fe_0 + tr_y_xy_xy[i] * pa_x[i];

        tr_y_xxy_xz[i] = ts_xx_xz[i] * fe_0 + tr_y_xx_xz[i] * pa_y[i];

        tr_y_xxy_yy[i] = tr_y_y_yy[i] * fe_0 + tr_y_xy_yy[i] * pa_x[i];

        tr_y_xxy_yz[i] = tr_y_y_yz[i] * fe_0 + tr_y_xy_yz[i] * pa_x[i];

        tr_y_xxy_zz[i] = tr_y_y_zz[i] * fe_0 + tr_y_xy_zz[i] * pa_x[i];
    }

    // Set up 72-78 components of targeted buffer : FD

    auto tr_y_xxz_xx = pbuffer.data(idx_dip_fd + 72);

    auto tr_y_xxz_xy = pbuffer.data(idx_dip_fd + 73);

    auto tr_y_xxz_xz = pbuffer.data(idx_dip_fd + 74);

    auto tr_y_xxz_yy = pbuffer.data(idx_dip_fd + 75);

    auto tr_y_xxz_yz = pbuffer.data(idx_dip_fd + 76);

    auto tr_y_xxz_zz = pbuffer.data(idx_dip_fd + 77);

#pragma omp simd aligned(pa_x,            \
                             pa_z,        \
                             tr_y_xx_x,   \
                             tr_y_xx_xx,  \
                             tr_y_xx_xy,  \
                             tr_y_xx_xz,  \
                             tr_y_xx_yy,  \
                             tr_y_xxz_xx, \
                             tr_y_xxz_xy, \
                             tr_y_xxz_xz, \
                             tr_y_xxz_yy, \
                             tr_y_xxz_yz, \
                             tr_y_xxz_zz, \
                             tr_y_xz_yz,  \
                             tr_y_xz_zz,  \
                             tr_y_z_yz,   \
                             tr_y_z_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxz_xx[i] = tr_y_xx_xx[i] * pa_z[i];

        tr_y_xxz_xy[i] = tr_y_xx_xy[i] * pa_z[i];

        tr_y_xxz_xz[i] = tr_y_xx_x[i] * fe_0 + tr_y_xx_xz[i] * pa_z[i];

        tr_y_xxz_yy[i] = tr_y_xx_yy[i] * pa_z[i];

        tr_y_xxz_yz[i] = tr_y_z_yz[i] * fe_0 + tr_y_xz_yz[i] * pa_x[i];

        tr_y_xxz_zz[i] = tr_y_z_zz[i] * fe_0 + tr_y_xz_zz[i] * pa_x[i];
    }

    // Set up 78-84 components of targeted buffer : FD

    auto tr_y_xyy_xx = pbuffer.data(idx_dip_fd + 78);

    auto tr_y_xyy_xy = pbuffer.data(idx_dip_fd + 79);

    auto tr_y_xyy_xz = pbuffer.data(idx_dip_fd + 80);

    auto tr_y_xyy_yy = pbuffer.data(idx_dip_fd + 81);

    auto tr_y_xyy_yz = pbuffer.data(idx_dip_fd + 82);

    auto tr_y_xyy_zz = pbuffer.data(idx_dip_fd + 83);

#pragma omp simd aligned(pa_x,            \
                             tr_y_xyy_xx, \
                             tr_y_xyy_xy, \
                             tr_y_xyy_xz, \
                             tr_y_xyy_yy, \
                             tr_y_xyy_yz, \
                             tr_y_xyy_zz, \
                             tr_y_yy_x,   \
                             tr_y_yy_xx,  \
                             tr_y_yy_xy,  \
                             tr_y_yy_xz,  \
                             tr_y_yy_y,   \
                             tr_y_yy_yy,  \
                             tr_y_yy_yz,  \
                             tr_y_yy_z,   \
                             tr_y_yy_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyy_xx[i] = 2.0 * tr_y_yy_x[i] * fe_0 + tr_y_yy_xx[i] * pa_x[i];

        tr_y_xyy_xy[i] = tr_y_yy_y[i] * fe_0 + tr_y_yy_xy[i] * pa_x[i];

        tr_y_xyy_xz[i] = tr_y_yy_z[i] * fe_0 + tr_y_yy_xz[i] * pa_x[i];

        tr_y_xyy_yy[i] = tr_y_yy_yy[i] * pa_x[i];

        tr_y_xyy_yz[i] = tr_y_yy_yz[i] * pa_x[i];

        tr_y_xyy_zz[i] = tr_y_yy_zz[i] * pa_x[i];
    }

    // Set up 84-90 components of targeted buffer : FD

    auto tr_y_xyz_xx = pbuffer.data(idx_dip_fd + 84);

    auto tr_y_xyz_xy = pbuffer.data(idx_dip_fd + 85);

    auto tr_y_xyz_xz = pbuffer.data(idx_dip_fd + 86);

    auto tr_y_xyz_yy = pbuffer.data(idx_dip_fd + 87);

    auto tr_y_xyz_yz = pbuffer.data(idx_dip_fd + 88);

    auto tr_y_xyz_zz = pbuffer.data(idx_dip_fd + 89);

#pragma omp simd aligned(pa_x,            \
                             pa_z,        \
                             tr_y_xy_xx,  \
                             tr_y_xy_xy,  \
                             tr_y_xyz_xx, \
                             tr_y_xyz_xy, \
                             tr_y_xyz_xz, \
                             tr_y_xyz_yy, \
                             tr_y_xyz_yz, \
                             tr_y_xyz_zz, \
                             tr_y_yz_xz,  \
                             tr_y_yz_yy,  \
                             tr_y_yz_yz,  \
                             tr_y_yz_z,   \
                             tr_y_yz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyz_xx[i] = tr_y_xy_xx[i] * pa_z[i];

        tr_y_xyz_xy[i] = tr_y_xy_xy[i] * pa_z[i];

        tr_y_xyz_xz[i] = tr_y_yz_z[i] * fe_0 + tr_y_yz_xz[i] * pa_x[i];

        tr_y_xyz_yy[i] = tr_y_yz_yy[i] * pa_x[i];

        tr_y_xyz_yz[i] = tr_y_yz_yz[i] * pa_x[i];

        tr_y_xyz_zz[i] = tr_y_yz_zz[i] * pa_x[i];
    }

    // Set up 90-96 components of targeted buffer : FD

    auto tr_y_xzz_xx = pbuffer.data(idx_dip_fd + 90);

    auto tr_y_xzz_xy = pbuffer.data(idx_dip_fd + 91);

    auto tr_y_xzz_xz = pbuffer.data(idx_dip_fd + 92);

    auto tr_y_xzz_yy = pbuffer.data(idx_dip_fd + 93);

    auto tr_y_xzz_yz = pbuffer.data(idx_dip_fd + 94);

    auto tr_y_xzz_zz = pbuffer.data(idx_dip_fd + 95);

#pragma omp simd aligned(pa_x,            \
                             tr_y_xzz_xx, \
                             tr_y_xzz_xy, \
                             tr_y_xzz_xz, \
                             tr_y_xzz_yy, \
                             tr_y_xzz_yz, \
                             tr_y_xzz_zz, \
                             tr_y_zz_x,   \
                             tr_y_zz_xx,  \
                             tr_y_zz_xy,  \
                             tr_y_zz_xz,  \
                             tr_y_zz_y,   \
                             tr_y_zz_yy,  \
                             tr_y_zz_yz,  \
                             tr_y_zz_z,   \
                             tr_y_zz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xzz_xx[i] = 2.0 * tr_y_zz_x[i] * fe_0 + tr_y_zz_xx[i] * pa_x[i];

        tr_y_xzz_xy[i] = tr_y_zz_y[i] * fe_0 + tr_y_zz_xy[i] * pa_x[i];

        tr_y_xzz_xz[i] = tr_y_zz_z[i] * fe_0 + tr_y_zz_xz[i] * pa_x[i];

        tr_y_xzz_yy[i] = tr_y_zz_yy[i] * pa_x[i];

        tr_y_xzz_yz[i] = tr_y_zz_yz[i] * pa_x[i];

        tr_y_xzz_zz[i] = tr_y_zz_zz[i] * pa_x[i];
    }

    // Set up 96-102 components of targeted buffer : FD

    auto tr_y_yyy_xx = pbuffer.data(idx_dip_fd + 96);

    auto tr_y_yyy_xy = pbuffer.data(idx_dip_fd + 97);

    auto tr_y_yyy_xz = pbuffer.data(idx_dip_fd + 98);

    auto tr_y_yyy_yy = pbuffer.data(idx_dip_fd + 99);

    auto tr_y_yyy_yz = pbuffer.data(idx_dip_fd + 100);

    auto tr_y_yyy_zz = pbuffer.data(idx_dip_fd + 101);

#pragma omp simd aligned(pa_y,            \
                             tr_y_y_xx,   \
                             tr_y_y_xy,   \
                             tr_y_y_xz,   \
                             tr_y_y_yy,   \
                             tr_y_y_yz,   \
                             tr_y_y_zz,   \
                             tr_y_yy_x,   \
                             tr_y_yy_xx,  \
                             tr_y_yy_xy,  \
                             tr_y_yy_xz,  \
                             tr_y_yy_y,   \
                             tr_y_yy_yy,  \
                             tr_y_yy_yz,  \
                             tr_y_yy_z,   \
                             tr_y_yy_zz,  \
                             tr_y_yyy_xx, \
                             tr_y_yyy_xy, \
                             tr_y_yyy_xz, \
                             tr_y_yyy_yy, \
                             tr_y_yyy_yz, \
                             tr_y_yyy_zz, \
                             ts_yy_xx,    \
                             ts_yy_xy,    \
                             ts_yy_xz,    \
                             ts_yy_yy,    \
                             ts_yy_yz,    \
                             ts_yy_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyy_xx[i] = 2.0 * tr_y_y_xx[i] * fe_0 + ts_yy_xx[i] * fe_0 + tr_y_yy_xx[i] * pa_y[i];

        tr_y_yyy_xy[i] = 2.0 * tr_y_y_xy[i] * fe_0 + tr_y_yy_x[i] * fe_0 + ts_yy_xy[i] * fe_0 + tr_y_yy_xy[i] * pa_y[i];

        tr_y_yyy_xz[i] = 2.0 * tr_y_y_xz[i] * fe_0 + ts_yy_xz[i] * fe_0 + tr_y_yy_xz[i] * pa_y[i];

        tr_y_yyy_yy[i] = 2.0 * tr_y_y_yy[i] * fe_0 + 2.0 * tr_y_yy_y[i] * fe_0 + ts_yy_yy[i] * fe_0 + tr_y_yy_yy[i] * pa_y[i];

        tr_y_yyy_yz[i] = 2.0 * tr_y_y_yz[i] * fe_0 + tr_y_yy_z[i] * fe_0 + ts_yy_yz[i] * fe_0 + tr_y_yy_yz[i] * pa_y[i];

        tr_y_yyy_zz[i] = 2.0 * tr_y_y_zz[i] * fe_0 + ts_yy_zz[i] * fe_0 + tr_y_yy_zz[i] * pa_y[i];
    }

    // Set up 102-108 components of targeted buffer : FD

    auto tr_y_yyz_xx = pbuffer.data(idx_dip_fd + 102);

    auto tr_y_yyz_xy = pbuffer.data(idx_dip_fd + 103);

    auto tr_y_yyz_xz = pbuffer.data(idx_dip_fd + 104);

    auto tr_y_yyz_yy = pbuffer.data(idx_dip_fd + 105);

    auto tr_y_yyz_yz = pbuffer.data(idx_dip_fd + 106);

    auto tr_y_yyz_zz = pbuffer.data(idx_dip_fd + 107);

#pragma omp simd aligned(pa_z,            \
                             tr_y_yy_x,   \
                             tr_y_yy_xx,  \
                             tr_y_yy_xy,  \
                             tr_y_yy_xz,  \
                             tr_y_yy_y,   \
                             tr_y_yy_yy,  \
                             tr_y_yy_yz,  \
                             tr_y_yy_z,   \
                             tr_y_yy_zz,  \
                             tr_y_yyz_xx, \
                             tr_y_yyz_xy, \
                             tr_y_yyz_xz, \
                             tr_y_yyz_yy, \
                             tr_y_yyz_yz, \
                             tr_y_yyz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyz_xx[i] = tr_y_yy_xx[i] * pa_z[i];

        tr_y_yyz_xy[i] = tr_y_yy_xy[i] * pa_z[i];

        tr_y_yyz_xz[i] = tr_y_yy_x[i] * fe_0 + tr_y_yy_xz[i] * pa_z[i];

        tr_y_yyz_yy[i] = tr_y_yy_yy[i] * pa_z[i];

        tr_y_yyz_yz[i] = tr_y_yy_y[i] * fe_0 + tr_y_yy_yz[i] * pa_z[i];

        tr_y_yyz_zz[i] = 2.0 * tr_y_yy_z[i] * fe_0 + tr_y_yy_zz[i] * pa_z[i];
    }

    // Set up 108-114 components of targeted buffer : FD

    auto tr_y_yzz_xx = pbuffer.data(idx_dip_fd + 108);

    auto tr_y_yzz_xy = pbuffer.data(idx_dip_fd + 109);

    auto tr_y_yzz_xz = pbuffer.data(idx_dip_fd + 110);

    auto tr_y_yzz_yy = pbuffer.data(idx_dip_fd + 111);

    auto tr_y_yzz_yz = pbuffer.data(idx_dip_fd + 112);

    auto tr_y_yzz_zz = pbuffer.data(idx_dip_fd + 113);

#pragma omp simd aligned(pa_y,            \
                             pa_z,        \
                             tr_y_y_xy,   \
                             tr_y_y_yy,   \
                             tr_y_yz_xy,  \
                             tr_y_yz_yy,  \
                             tr_y_yzz_xx, \
                             tr_y_yzz_xy, \
                             tr_y_yzz_xz, \
                             tr_y_yzz_yy, \
                             tr_y_yzz_yz, \
                             tr_y_yzz_zz, \
                             tr_y_zz_xx,  \
                             tr_y_zz_xz,  \
                             tr_y_zz_yz,  \
                             tr_y_zz_z,   \
                             tr_y_zz_zz,  \
                             ts_zz_xx,    \
                             ts_zz_xz,    \
                             ts_zz_yz,    \
                             ts_zz_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yzz_xx[i] = ts_zz_xx[i] * fe_0 + tr_y_zz_xx[i] * pa_y[i];

        tr_y_yzz_xy[i] = tr_y_y_xy[i] * fe_0 + tr_y_yz_xy[i] * pa_z[i];

        tr_y_yzz_xz[i] = ts_zz_xz[i] * fe_0 + tr_y_zz_xz[i] * pa_y[i];

        tr_y_yzz_yy[i] = tr_y_y_yy[i] * fe_0 + tr_y_yz_yy[i] * pa_z[i];

        tr_y_yzz_yz[i] = tr_y_zz_z[i] * fe_0 + ts_zz_yz[i] * fe_0 + tr_y_zz_yz[i] * pa_y[i];

        tr_y_yzz_zz[i] = ts_zz_zz[i] * fe_0 + tr_y_zz_zz[i] * pa_y[i];
    }

    // Set up 114-120 components of targeted buffer : FD

    auto tr_y_zzz_xx = pbuffer.data(idx_dip_fd + 114);

    auto tr_y_zzz_xy = pbuffer.data(idx_dip_fd + 115);

    auto tr_y_zzz_xz = pbuffer.data(idx_dip_fd + 116);

    auto tr_y_zzz_yy = pbuffer.data(idx_dip_fd + 117);

    auto tr_y_zzz_yz = pbuffer.data(idx_dip_fd + 118);

    auto tr_y_zzz_zz = pbuffer.data(idx_dip_fd + 119);

#pragma omp simd aligned(pa_z,            \
                             tr_y_z_xx,   \
                             tr_y_z_xy,   \
                             tr_y_z_xz,   \
                             tr_y_z_yy,   \
                             tr_y_z_yz,   \
                             tr_y_z_zz,   \
                             tr_y_zz_x,   \
                             tr_y_zz_xx,  \
                             tr_y_zz_xy,  \
                             tr_y_zz_xz,  \
                             tr_y_zz_y,   \
                             tr_y_zz_yy,  \
                             tr_y_zz_yz,  \
                             tr_y_zz_z,   \
                             tr_y_zz_zz,  \
                             tr_y_zzz_xx, \
                             tr_y_zzz_xy, \
                             tr_y_zzz_xz, \
                             tr_y_zzz_yy, \
                             tr_y_zzz_yz, \
                             tr_y_zzz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_zzz_xx[i] = 2.0 * tr_y_z_xx[i] * fe_0 + tr_y_zz_xx[i] * pa_z[i];

        tr_y_zzz_xy[i] = 2.0 * tr_y_z_xy[i] * fe_0 + tr_y_zz_xy[i] * pa_z[i];

        tr_y_zzz_xz[i] = 2.0 * tr_y_z_xz[i] * fe_0 + tr_y_zz_x[i] * fe_0 + tr_y_zz_xz[i] * pa_z[i];

        tr_y_zzz_yy[i] = 2.0 * tr_y_z_yy[i] * fe_0 + tr_y_zz_yy[i] * pa_z[i];

        tr_y_zzz_yz[i] = 2.0 * tr_y_z_yz[i] * fe_0 + tr_y_zz_y[i] * fe_0 + tr_y_zz_yz[i] * pa_z[i];

        tr_y_zzz_zz[i] = 2.0 * tr_y_z_zz[i] * fe_0 + 2.0 * tr_y_zz_z[i] * fe_0 + tr_y_zz_zz[i] * pa_z[i];
    }

    // Set up 120-126 components of targeted buffer : FD

    auto tr_z_xxx_xx = pbuffer.data(idx_dip_fd + 120);

    auto tr_z_xxx_xy = pbuffer.data(idx_dip_fd + 121);

    auto tr_z_xxx_xz = pbuffer.data(idx_dip_fd + 122);

    auto tr_z_xxx_yy = pbuffer.data(idx_dip_fd + 123);

    auto tr_z_xxx_yz = pbuffer.data(idx_dip_fd + 124);

    auto tr_z_xxx_zz = pbuffer.data(idx_dip_fd + 125);

#pragma omp simd aligned(pa_x,            \
                             tr_z_x_xx,   \
                             tr_z_x_xy,   \
                             tr_z_x_xz,   \
                             tr_z_x_yy,   \
                             tr_z_x_yz,   \
                             tr_z_x_zz,   \
                             tr_z_xx_x,   \
                             tr_z_xx_xx,  \
                             tr_z_xx_xy,  \
                             tr_z_xx_xz,  \
                             tr_z_xx_y,   \
                             tr_z_xx_yy,  \
                             tr_z_xx_yz,  \
                             tr_z_xx_z,   \
                             tr_z_xx_zz,  \
                             tr_z_xxx_xx, \
                             tr_z_xxx_xy, \
                             tr_z_xxx_xz, \
                             tr_z_xxx_yy, \
                             tr_z_xxx_yz, \
                             tr_z_xxx_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxx_xx[i] = 2.0 * tr_z_x_xx[i] * fe_0 + 2.0 * tr_z_xx_x[i] * fe_0 + tr_z_xx_xx[i] * pa_x[i];

        tr_z_xxx_xy[i] = 2.0 * tr_z_x_xy[i] * fe_0 + tr_z_xx_y[i] * fe_0 + tr_z_xx_xy[i] * pa_x[i];

        tr_z_xxx_xz[i] = 2.0 * tr_z_x_xz[i] * fe_0 + tr_z_xx_z[i] * fe_0 + tr_z_xx_xz[i] * pa_x[i];

        tr_z_xxx_yy[i] = 2.0 * tr_z_x_yy[i] * fe_0 + tr_z_xx_yy[i] * pa_x[i];

        tr_z_xxx_yz[i] = 2.0 * tr_z_x_yz[i] * fe_0 + tr_z_xx_yz[i] * pa_x[i];

        tr_z_xxx_zz[i] = 2.0 * tr_z_x_zz[i] * fe_0 + tr_z_xx_zz[i] * pa_x[i];
    }

    // Set up 126-132 components of targeted buffer : FD

    auto tr_z_xxy_xx = pbuffer.data(idx_dip_fd + 126);

    auto tr_z_xxy_xy = pbuffer.data(idx_dip_fd + 127);

    auto tr_z_xxy_xz = pbuffer.data(idx_dip_fd + 128);

    auto tr_z_xxy_yy = pbuffer.data(idx_dip_fd + 129);

    auto tr_z_xxy_yz = pbuffer.data(idx_dip_fd + 130);

    auto tr_z_xxy_zz = pbuffer.data(idx_dip_fd + 131);

#pragma omp simd aligned(pa_x,            \
                             pa_y,        \
                             tr_z_xx_x,   \
                             tr_z_xx_xx,  \
                             tr_z_xx_xy,  \
                             tr_z_xx_xz,  \
                             tr_z_xx_zz,  \
                             tr_z_xxy_xx, \
                             tr_z_xxy_xy, \
                             tr_z_xxy_xz, \
                             tr_z_xxy_yy, \
                             tr_z_xxy_yz, \
                             tr_z_xxy_zz, \
                             tr_z_xy_yy,  \
                             tr_z_xy_yz,  \
                             tr_z_y_yy,   \
                             tr_z_y_yz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxy_xx[i] = tr_z_xx_xx[i] * pa_y[i];

        tr_z_xxy_xy[i] = tr_z_xx_x[i] * fe_0 + tr_z_xx_xy[i] * pa_y[i];

        tr_z_xxy_xz[i] = tr_z_xx_xz[i] * pa_y[i];

        tr_z_xxy_yy[i] = tr_z_y_yy[i] * fe_0 + tr_z_xy_yy[i] * pa_x[i];

        tr_z_xxy_yz[i] = tr_z_y_yz[i] * fe_0 + tr_z_xy_yz[i] * pa_x[i];

        tr_z_xxy_zz[i] = tr_z_xx_zz[i] * pa_y[i];
    }

    // Set up 132-138 components of targeted buffer : FD

    auto tr_z_xxz_xx = pbuffer.data(idx_dip_fd + 132);

    auto tr_z_xxz_xy = pbuffer.data(idx_dip_fd + 133);

    auto tr_z_xxz_xz = pbuffer.data(idx_dip_fd + 134);

    auto tr_z_xxz_yy = pbuffer.data(idx_dip_fd + 135);

    auto tr_z_xxz_yz = pbuffer.data(idx_dip_fd + 136);

    auto tr_z_xxz_zz = pbuffer.data(idx_dip_fd + 137);

#pragma omp simd aligned(pa_x,            \
                             pa_z,        \
                             tr_z_xx_xx,  \
                             tr_z_xx_xy,  \
                             tr_z_xxz_xx, \
                             tr_z_xxz_xy, \
                             tr_z_xxz_xz, \
                             tr_z_xxz_yy, \
                             tr_z_xxz_yz, \
                             tr_z_xxz_zz, \
                             tr_z_xz_xz,  \
                             tr_z_xz_yy,  \
                             tr_z_xz_yz,  \
                             tr_z_xz_z,   \
                             tr_z_xz_zz,  \
                             tr_z_z_xz,   \
                             tr_z_z_yy,   \
                             tr_z_z_yz,   \
                             tr_z_z_zz,   \
                             ts_xx_xx,    \
                             ts_xx_xy,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxz_xx[i] = ts_xx_xx[i] * fe_0 + tr_z_xx_xx[i] * pa_z[i];

        tr_z_xxz_xy[i] = ts_xx_xy[i] * fe_0 + tr_z_xx_xy[i] * pa_z[i];

        tr_z_xxz_xz[i] = tr_z_z_xz[i] * fe_0 + tr_z_xz_z[i] * fe_0 + tr_z_xz_xz[i] * pa_x[i];

        tr_z_xxz_yy[i] = tr_z_z_yy[i] * fe_0 + tr_z_xz_yy[i] * pa_x[i];

        tr_z_xxz_yz[i] = tr_z_z_yz[i] * fe_0 + tr_z_xz_yz[i] * pa_x[i];

        tr_z_xxz_zz[i] = tr_z_z_zz[i] * fe_0 + tr_z_xz_zz[i] * pa_x[i];
    }

    // Set up 138-144 components of targeted buffer : FD

    auto tr_z_xyy_xx = pbuffer.data(idx_dip_fd + 138);

    auto tr_z_xyy_xy = pbuffer.data(idx_dip_fd + 139);

    auto tr_z_xyy_xz = pbuffer.data(idx_dip_fd + 140);

    auto tr_z_xyy_yy = pbuffer.data(idx_dip_fd + 141);

    auto tr_z_xyy_yz = pbuffer.data(idx_dip_fd + 142);

    auto tr_z_xyy_zz = pbuffer.data(idx_dip_fd + 143);

#pragma omp simd aligned(pa_x,            \
                             tr_z_xyy_xx, \
                             tr_z_xyy_xy, \
                             tr_z_xyy_xz, \
                             tr_z_xyy_yy, \
                             tr_z_xyy_yz, \
                             tr_z_xyy_zz, \
                             tr_z_yy_x,   \
                             tr_z_yy_xx,  \
                             tr_z_yy_xy,  \
                             tr_z_yy_xz,  \
                             tr_z_yy_y,   \
                             tr_z_yy_yy,  \
                             tr_z_yy_yz,  \
                             tr_z_yy_z,   \
                             tr_z_yy_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyy_xx[i] = 2.0 * tr_z_yy_x[i] * fe_0 + tr_z_yy_xx[i] * pa_x[i];

        tr_z_xyy_xy[i] = tr_z_yy_y[i] * fe_0 + tr_z_yy_xy[i] * pa_x[i];

        tr_z_xyy_xz[i] = tr_z_yy_z[i] * fe_0 + tr_z_yy_xz[i] * pa_x[i];

        tr_z_xyy_yy[i] = tr_z_yy_yy[i] * pa_x[i];

        tr_z_xyy_yz[i] = tr_z_yy_yz[i] * pa_x[i];

        tr_z_xyy_zz[i] = tr_z_yy_zz[i] * pa_x[i];
    }

    // Set up 144-150 components of targeted buffer : FD

    auto tr_z_xyz_xx = pbuffer.data(idx_dip_fd + 144);

    auto tr_z_xyz_xy = pbuffer.data(idx_dip_fd + 145);

    auto tr_z_xyz_xz = pbuffer.data(idx_dip_fd + 146);

    auto tr_z_xyz_yy = pbuffer.data(idx_dip_fd + 147);

    auto tr_z_xyz_yz = pbuffer.data(idx_dip_fd + 148);

    auto tr_z_xyz_zz = pbuffer.data(idx_dip_fd + 149);

#pragma omp simd aligned(pa_x,            \
                             pa_y,        \
                             tr_z_xyz_xx, \
                             tr_z_xyz_xy, \
                             tr_z_xyz_xz, \
                             tr_z_xyz_yy, \
                             tr_z_xyz_yz, \
                             tr_z_xyz_zz, \
                             tr_z_xz_xx,  \
                             tr_z_xz_xz,  \
                             tr_z_yz_xy,  \
                             tr_z_yz_y,   \
                             tr_z_yz_yy,  \
                             tr_z_yz_yz,  \
                             tr_z_yz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyz_xx[i] = tr_z_xz_xx[i] * pa_y[i];

        tr_z_xyz_xy[i] = tr_z_yz_y[i] * fe_0 + tr_z_yz_xy[i] * pa_x[i];

        tr_z_xyz_xz[i] = tr_z_xz_xz[i] * pa_y[i];

        tr_z_xyz_yy[i] = tr_z_yz_yy[i] * pa_x[i];

        tr_z_xyz_yz[i] = tr_z_yz_yz[i] * pa_x[i];

        tr_z_xyz_zz[i] = tr_z_yz_zz[i] * pa_x[i];
    }

    // Set up 150-156 components of targeted buffer : FD

    auto tr_z_xzz_xx = pbuffer.data(idx_dip_fd + 150);

    auto tr_z_xzz_xy = pbuffer.data(idx_dip_fd + 151);

    auto tr_z_xzz_xz = pbuffer.data(idx_dip_fd + 152);

    auto tr_z_xzz_yy = pbuffer.data(idx_dip_fd + 153);

    auto tr_z_xzz_yz = pbuffer.data(idx_dip_fd + 154);

    auto tr_z_xzz_zz = pbuffer.data(idx_dip_fd + 155);

#pragma omp simd aligned(pa_x,            \
                             tr_z_xzz_xx, \
                             tr_z_xzz_xy, \
                             tr_z_xzz_xz, \
                             tr_z_xzz_yy, \
                             tr_z_xzz_yz, \
                             tr_z_xzz_zz, \
                             tr_z_zz_x,   \
                             tr_z_zz_xx,  \
                             tr_z_zz_xy,  \
                             tr_z_zz_xz,  \
                             tr_z_zz_y,   \
                             tr_z_zz_yy,  \
                             tr_z_zz_yz,  \
                             tr_z_zz_z,   \
                             tr_z_zz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xzz_xx[i] = 2.0 * tr_z_zz_x[i] * fe_0 + tr_z_zz_xx[i] * pa_x[i];

        tr_z_xzz_xy[i] = tr_z_zz_y[i] * fe_0 + tr_z_zz_xy[i] * pa_x[i];

        tr_z_xzz_xz[i] = tr_z_zz_z[i] * fe_0 + tr_z_zz_xz[i] * pa_x[i];

        tr_z_xzz_yy[i] = tr_z_zz_yy[i] * pa_x[i];

        tr_z_xzz_yz[i] = tr_z_zz_yz[i] * pa_x[i];

        tr_z_xzz_zz[i] = tr_z_zz_zz[i] * pa_x[i];
    }

    // Set up 156-162 components of targeted buffer : FD

    auto tr_z_yyy_xx = pbuffer.data(idx_dip_fd + 156);

    auto tr_z_yyy_xy = pbuffer.data(idx_dip_fd + 157);

    auto tr_z_yyy_xz = pbuffer.data(idx_dip_fd + 158);

    auto tr_z_yyy_yy = pbuffer.data(idx_dip_fd + 159);

    auto tr_z_yyy_yz = pbuffer.data(idx_dip_fd + 160);

    auto tr_z_yyy_zz = pbuffer.data(idx_dip_fd + 161);

#pragma omp simd aligned(pa_y,            \
                             tr_z_y_xx,   \
                             tr_z_y_xy,   \
                             tr_z_y_xz,   \
                             tr_z_y_yy,   \
                             tr_z_y_yz,   \
                             tr_z_y_zz,   \
                             tr_z_yy_x,   \
                             tr_z_yy_xx,  \
                             tr_z_yy_xy,  \
                             tr_z_yy_xz,  \
                             tr_z_yy_y,   \
                             tr_z_yy_yy,  \
                             tr_z_yy_yz,  \
                             tr_z_yy_z,   \
                             tr_z_yy_zz,  \
                             tr_z_yyy_xx, \
                             tr_z_yyy_xy, \
                             tr_z_yyy_xz, \
                             tr_z_yyy_yy, \
                             tr_z_yyy_yz, \
                             tr_z_yyy_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyy_xx[i] = 2.0 * tr_z_y_xx[i] * fe_0 + tr_z_yy_xx[i] * pa_y[i];

        tr_z_yyy_xy[i] = 2.0 * tr_z_y_xy[i] * fe_0 + tr_z_yy_x[i] * fe_0 + tr_z_yy_xy[i] * pa_y[i];

        tr_z_yyy_xz[i] = 2.0 * tr_z_y_xz[i] * fe_0 + tr_z_yy_xz[i] * pa_y[i];

        tr_z_yyy_yy[i] = 2.0 * tr_z_y_yy[i] * fe_0 + 2.0 * tr_z_yy_y[i] * fe_0 + tr_z_yy_yy[i] * pa_y[i];

        tr_z_yyy_yz[i] = 2.0 * tr_z_y_yz[i] * fe_0 + tr_z_yy_z[i] * fe_0 + tr_z_yy_yz[i] * pa_y[i];

        tr_z_yyy_zz[i] = 2.0 * tr_z_y_zz[i] * fe_0 + tr_z_yy_zz[i] * pa_y[i];
    }

    // Set up 162-168 components of targeted buffer : FD

    auto tr_z_yyz_xx = pbuffer.data(idx_dip_fd + 162);

    auto tr_z_yyz_xy = pbuffer.data(idx_dip_fd + 163);

    auto tr_z_yyz_xz = pbuffer.data(idx_dip_fd + 164);

    auto tr_z_yyz_yy = pbuffer.data(idx_dip_fd + 165);

    auto tr_z_yyz_yz = pbuffer.data(idx_dip_fd + 166);

    auto tr_z_yyz_zz = pbuffer.data(idx_dip_fd + 167);

#pragma omp simd aligned(pa_y,            \
                             pa_z,        \
                             tr_z_yy_xy,  \
                             tr_z_yy_yy,  \
                             tr_z_yyz_xx, \
                             tr_z_yyz_xy, \
                             tr_z_yyz_xz, \
                             tr_z_yyz_yy, \
                             tr_z_yyz_yz, \
                             tr_z_yyz_zz, \
                             tr_z_yz_xx,  \
                             tr_z_yz_xz,  \
                             tr_z_yz_yz,  \
                             tr_z_yz_z,   \
                             tr_z_yz_zz,  \
                             tr_z_z_xx,   \
                             tr_z_z_xz,   \
                             tr_z_z_yz,   \
                             tr_z_z_zz,   \
                             ts_yy_xy,    \
                             ts_yy_yy,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyz_xx[i] = tr_z_z_xx[i] * fe_0 + tr_z_yz_xx[i] * pa_y[i];

        tr_z_yyz_xy[i] = ts_yy_xy[i] * fe_0 + tr_z_yy_xy[i] * pa_z[i];

        tr_z_yyz_xz[i] = tr_z_z_xz[i] * fe_0 + tr_z_yz_xz[i] * pa_y[i];

        tr_z_yyz_yy[i] = ts_yy_yy[i] * fe_0 + tr_z_yy_yy[i] * pa_z[i];

        tr_z_yyz_yz[i] = tr_z_z_yz[i] * fe_0 + tr_z_yz_z[i] * fe_0 + tr_z_yz_yz[i] * pa_y[i];

        tr_z_yyz_zz[i] = tr_z_z_zz[i] * fe_0 + tr_z_yz_zz[i] * pa_y[i];
    }

    // Set up 168-174 components of targeted buffer : FD

    auto tr_z_yzz_xx = pbuffer.data(idx_dip_fd + 168);

    auto tr_z_yzz_xy = pbuffer.data(idx_dip_fd + 169);

    auto tr_z_yzz_xz = pbuffer.data(idx_dip_fd + 170);

    auto tr_z_yzz_yy = pbuffer.data(idx_dip_fd + 171);

    auto tr_z_yzz_yz = pbuffer.data(idx_dip_fd + 172);

    auto tr_z_yzz_zz = pbuffer.data(idx_dip_fd + 173);

#pragma omp simd aligned(pa_y,            \
                             tr_z_yzz_xx, \
                             tr_z_yzz_xy, \
                             tr_z_yzz_xz, \
                             tr_z_yzz_yy, \
                             tr_z_yzz_yz, \
                             tr_z_yzz_zz, \
                             tr_z_zz_x,   \
                             tr_z_zz_xx,  \
                             tr_z_zz_xy,  \
                             tr_z_zz_xz,  \
                             tr_z_zz_y,   \
                             tr_z_zz_yy,  \
                             tr_z_zz_yz,  \
                             tr_z_zz_z,   \
                             tr_z_zz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yzz_xx[i] = tr_z_zz_xx[i] * pa_y[i];

        tr_z_yzz_xy[i] = tr_z_zz_x[i] * fe_0 + tr_z_zz_xy[i] * pa_y[i];

        tr_z_yzz_xz[i] = tr_z_zz_xz[i] * pa_y[i];

        tr_z_yzz_yy[i] = 2.0 * tr_z_zz_y[i] * fe_0 + tr_z_zz_yy[i] * pa_y[i];

        tr_z_yzz_yz[i] = tr_z_zz_z[i] * fe_0 + tr_z_zz_yz[i] * pa_y[i];

        tr_z_yzz_zz[i] = tr_z_zz_zz[i] * pa_y[i];
    }

    // Set up 174-180 components of targeted buffer : FD

    auto tr_z_zzz_xx = pbuffer.data(idx_dip_fd + 174);

    auto tr_z_zzz_xy = pbuffer.data(idx_dip_fd + 175);

    auto tr_z_zzz_xz = pbuffer.data(idx_dip_fd + 176);

    auto tr_z_zzz_yy = pbuffer.data(idx_dip_fd + 177);

    auto tr_z_zzz_yz = pbuffer.data(idx_dip_fd + 178);

    auto tr_z_zzz_zz = pbuffer.data(idx_dip_fd + 179);

#pragma omp simd aligned(pa_z,            \
                             tr_z_z_xx,   \
                             tr_z_z_xy,   \
                             tr_z_z_xz,   \
                             tr_z_z_yy,   \
                             tr_z_z_yz,   \
                             tr_z_z_zz,   \
                             tr_z_zz_x,   \
                             tr_z_zz_xx,  \
                             tr_z_zz_xy,  \
                             tr_z_zz_xz,  \
                             tr_z_zz_y,   \
                             tr_z_zz_yy,  \
                             tr_z_zz_yz,  \
                             tr_z_zz_z,   \
                             tr_z_zz_zz,  \
                             tr_z_zzz_xx, \
                             tr_z_zzz_xy, \
                             tr_z_zzz_xz, \
                             tr_z_zzz_yy, \
                             tr_z_zzz_yz, \
                             tr_z_zzz_zz, \
                             ts_zz_xx,    \
                             ts_zz_xy,    \
                             ts_zz_xz,    \
                             ts_zz_yy,    \
                             ts_zz_yz,    \
                             ts_zz_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_zzz_xx[i] = 2.0 * tr_z_z_xx[i] * fe_0 + ts_zz_xx[i] * fe_0 + tr_z_zz_xx[i] * pa_z[i];

        tr_z_zzz_xy[i] = 2.0 * tr_z_z_xy[i] * fe_0 + ts_zz_xy[i] * fe_0 + tr_z_zz_xy[i] * pa_z[i];

        tr_z_zzz_xz[i] = 2.0 * tr_z_z_xz[i] * fe_0 + tr_z_zz_x[i] * fe_0 + ts_zz_xz[i] * fe_0 + tr_z_zz_xz[i] * pa_z[i];

        tr_z_zzz_yy[i] = 2.0 * tr_z_z_yy[i] * fe_0 + ts_zz_yy[i] * fe_0 + tr_z_zz_yy[i] * pa_z[i];

        tr_z_zzz_yz[i] = 2.0 * tr_z_z_yz[i] * fe_0 + tr_z_zz_y[i] * fe_0 + ts_zz_yz[i] * fe_0 + tr_z_zz_yz[i] * pa_z[i];

        tr_z_zzz_zz[i] = 2.0 * tr_z_z_zz[i] * fe_0 + 2.0 * tr_z_zz_z[i] * fe_0 + ts_zz_zz[i] * fe_0 + tr_z_zz_zz[i] * pa_z[i];
    }
}

}  // namespace diprec
