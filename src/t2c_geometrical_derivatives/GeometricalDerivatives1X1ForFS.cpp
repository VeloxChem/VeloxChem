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

#include "GeometricalDerivatives1X1ForFS.hpp"

namespace t2cgeom {  // t2cgeom namespace

auto
comp_prim_op_geom_11_fs(CSimdArray<double>&       pbuffer,
                        const size_t              idx_op_geom_101_fs,
                        const size_t              idx_op_dp,
                        const size_t              idx_op_gp,
                        const size_t              op_comps,
                        const CSimdArray<double>& factors,
                        const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : DP

        auto to_xx_x = pbuffer.data(idx_op_dp + i * 18 + 0);

        auto to_xx_y = pbuffer.data(idx_op_dp + i * 18 + 1);

        auto to_xx_z = pbuffer.data(idx_op_dp + i * 18 + 2);

        auto to_xy_x = pbuffer.data(idx_op_dp + i * 18 + 3);

        auto to_xy_y = pbuffer.data(idx_op_dp + i * 18 + 4);

        auto to_xy_z = pbuffer.data(idx_op_dp + i * 18 + 5);

        auto to_xz_x = pbuffer.data(idx_op_dp + i * 18 + 6);

        auto to_xz_y = pbuffer.data(idx_op_dp + i * 18 + 7);

        auto to_xz_z = pbuffer.data(idx_op_dp + i * 18 + 8);

        auto to_yy_x = pbuffer.data(idx_op_dp + i * 18 + 9);

        auto to_yy_y = pbuffer.data(idx_op_dp + i * 18 + 10);

        auto to_yy_z = pbuffer.data(idx_op_dp + i * 18 + 11);

        auto to_yz_x = pbuffer.data(idx_op_dp + i * 18 + 12);

        auto to_yz_y = pbuffer.data(idx_op_dp + i * 18 + 13);

        auto to_yz_z = pbuffer.data(idx_op_dp + i * 18 + 14);

        auto to_zz_x = pbuffer.data(idx_op_dp + i * 18 + 15);

        auto to_zz_y = pbuffer.data(idx_op_dp + i * 18 + 16);

        auto to_zz_z = pbuffer.data(idx_op_dp + i * 18 + 17);

        // Set up components of auxiliary buffer : GP

        auto to_xxxx_x = pbuffer.data(idx_op_gp + i * 45 + 0);

        auto to_xxxx_y = pbuffer.data(idx_op_gp + i * 45 + 1);

        auto to_xxxx_z = pbuffer.data(idx_op_gp + i * 45 + 2);

        auto to_xxxy_x = pbuffer.data(idx_op_gp + i * 45 + 3);

        auto to_xxxy_y = pbuffer.data(idx_op_gp + i * 45 + 4);

        auto to_xxxy_z = pbuffer.data(idx_op_gp + i * 45 + 5);

        auto to_xxxz_x = pbuffer.data(idx_op_gp + i * 45 + 6);

        auto to_xxxz_y = pbuffer.data(idx_op_gp + i * 45 + 7);

        auto to_xxxz_z = pbuffer.data(idx_op_gp + i * 45 + 8);

        auto to_xxyy_x = pbuffer.data(idx_op_gp + i * 45 + 9);

        auto to_xxyy_y = pbuffer.data(idx_op_gp + i * 45 + 10);

        auto to_xxyy_z = pbuffer.data(idx_op_gp + i * 45 + 11);

        auto to_xxyz_x = pbuffer.data(idx_op_gp + i * 45 + 12);

        auto to_xxyz_y = pbuffer.data(idx_op_gp + i * 45 + 13);

        auto to_xxyz_z = pbuffer.data(idx_op_gp + i * 45 + 14);

        auto to_xxzz_x = pbuffer.data(idx_op_gp + i * 45 + 15);

        auto to_xxzz_y = pbuffer.data(idx_op_gp + i * 45 + 16);

        auto to_xxzz_z = pbuffer.data(idx_op_gp + i * 45 + 17);

        auto to_xyyy_x = pbuffer.data(idx_op_gp + i * 45 + 18);

        auto to_xyyy_y = pbuffer.data(idx_op_gp + i * 45 + 19);

        auto to_xyyy_z = pbuffer.data(idx_op_gp + i * 45 + 20);

        auto to_xyyz_x = pbuffer.data(idx_op_gp + i * 45 + 21);

        auto to_xyyz_y = pbuffer.data(idx_op_gp + i * 45 + 22);

        auto to_xyyz_z = pbuffer.data(idx_op_gp + i * 45 + 23);

        auto to_xyzz_x = pbuffer.data(idx_op_gp + i * 45 + 24);

        auto to_xyzz_y = pbuffer.data(idx_op_gp + i * 45 + 25);

        auto to_xyzz_z = pbuffer.data(idx_op_gp + i * 45 + 26);

        auto to_xzzz_x = pbuffer.data(idx_op_gp + i * 45 + 27);

        auto to_xzzz_y = pbuffer.data(idx_op_gp + i * 45 + 28);

        auto to_xzzz_z = pbuffer.data(idx_op_gp + i * 45 + 29);

        auto to_yyyy_x = pbuffer.data(idx_op_gp + i * 45 + 30);

        auto to_yyyy_y = pbuffer.data(idx_op_gp + i * 45 + 31);

        auto to_yyyy_z = pbuffer.data(idx_op_gp + i * 45 + 32);

        auto to_yyyz_x = pbuffer.data(idx_op_gp + i * 45 + 33);

        auto to_yyyz_y = pbuffer.data(idx_op_gp + i * 45 + 34);

        auto to_yyyz_z = pbuffer.data(idx_op_gp + i * 45 + 35);

        auto to_yyzz_x = pbuffer.data(idx_op_gp + i * 45 + 36);

        auto to_yyzz_y = pbuffer.data(idx_op_gp + i * 45 + 37);

        auto to_yyzz_z = pbuffer.data(idx_op_gp + i * 45 + 38);

        auto to_yzzz_x = pbuffer.data(idx_op_gp + i * 45 + 39);

        auto to_yzzz_y = pbuffer.data(idx_op_gp + i * 45 + 40);

        auto to_yzzz_z = pbuffer.data(idx_op_gp + i * 45 + 41);

        auto to_zzzz_x = pbuffer.data(idx_op_gp + i * 45 + 42);

        auto to_zzzz_y = pbuffer.data(idx_op_gp + i * 45 + 43);

        auto to_zzzz_z = pbuffer.data(idx_op_gp + i * 45 + 44);

        // Set up 0-10 components of targeted buffer : FS

        auto to_x_x_xxx_0 = pbuffer.data(idx_op_geom_101_fs + 0 * op_comps * 10 + i * 10 + 0);

        auto to_x_x_xxy_0 = pbuffer.data(idx_op_geom_101_fs + 0 * op_comps * 10 + i * 10 + 1);

        auto to_x_x_xxz_0 = pbuffer.data(idx_op_geom_101_fs + 0 * op_comps * 10 + i * 10 + 2);

        auto to_x_x_xyy_0 = pbuffer.data(idx_op_geom_101_fs + 0 * op_comps * 10 + i * 10 + 3);

        auto to_x_x_xyz_0 = pbuffer.data(idx_op_geom_101_fs + 0 * op_comps * 10 + i * 10 + 4);

        auto to_x_x_xzz_0 = pbuffer.data(idx_op_geom_101_fs + 0 * op_comps * 10 + i * 10 + 5);

        auto to_x_x_yyy_0 = pbuffer.data(idx_op_geom_101_fs + 0 * op_comps * 10 + i * 10 + 6);

        auto to_x_x_yyz_0 = pbuffer.data(idx_op_geom_101_fs + 0 * op_comps * 10 + i * 10 + 7);

        auto to_x_x_yzz_0 = pbuffer.data(idx_op_geom_101_fs + 0 * op_comps * 10 + i * 10 + 8);

        auto to_x_x_zzz_0 = pbuffer.data(idx_op_geom_101_fs + 0 * op_comps * 10 + i * 10 + 9);

#pragma omp simd aligned(to_x_x_xxx_0,     \
                             to_x_x_xxy_0, \
                             to_x_x_xxz_0, \
                             to_x_x_xyy_0, \
                             to_x_x_xyz_0, \
                             to_x_x_xzz_0, \
                             to_x_x_yyy_0, \
                             to_x_x_yyz_0, \
                             to_x_x_yzz_0, \
                             to_x_x_zzz_0, \
                             to_xx_x,      \
                             to_xxxx_x,    \
                             to_xxxy_x,    \
                             to_xxxz_x,    \
                             to_xxyy_x,    \
                             to_xxyz_x,    \
                             to_xxzz_x,    \
                             to_xy_x,      \
                             to_xyyy_x,    \
                             to_xyyz_x,    \
                             to_xyzz_x,    \
                             to_xz_x,      \
                             to_xzzz_x,    \
                             to_yy_x,      \
                             to_yz_x,      \
                             to_zz_x,      \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxx_0[k] = -6.0 * to_xx_x[k] * tke_0 + 4.0 * to_xxxx_x[k] * tbe_0 * tke_0;

            to_x_x_xxy_0[k] = -4.0 * to_xy_x[k] * tke_0 + 4.0 * to_xxxy_x[k] * tbe_0 * tke_0;

            to_x_x_xxz_0[k] = -4.0 * to_xz_x[k] * tke_0 + 4.0 * to_xxxz_x[k] * tbe_0 * tke_0;

            to_x_x_xyy_0[k] = -2.0 * to_yy_x[k] * tke_0 + 4.0 * to_xxyy_x[k] * tbe_0 * tke_0;

            to_x_x_xyz_0[k] = -2.0 * to_yz_x[k] * tke_0 + 4.0 * to_xxyz_x[k] * tbe_0 * tke_0;

            to_x_x_xzz_0[k] = -2.0 * to_zz_x[k] * tke_0 + 4.0 * to_xxzz_x[k] * tbe_0 * tke_0;

            to_x_x_yyy_0[k] = 4.0 * to_xyyy_x[k] * tbe_0 * tke_0;

            to_x_x_yyz_0[k] = 4.0 * to_xyyz_x[k] * tbe_0 * tke_0;

            to_x_x_yzz_0[k] = 4.0 * to_xyzz_x[k] * tbe_0 * tke_0;

            to_x_x_zzz_0[k] = 4.0 * to_xzzz_x[k] * tbe_0 * tke_0;
        }

        // Set up 10-20 components of targeted buffer : FS

        auto to_x_y_xxx_0 = pbuffer.data(idx_op_geom_101_fs + 1 * op_comps * 10 + i * 10 + 0);

        auto to_x_y_xxy_0 = pbuffer.data(idx_op_geom_101_fs + 1 * op_comps * 10 + i * 10 + 1);

        auto to_x_y_xxz_0 = pbuffer.data(idx_op_geom_101_fs + 1 * op_comps * 10 + i * 10 + 2);

        auto to_x_y_xyy_0 = pbuffer.data(idx_op_geom_101_fs + 1 * op_comps * 10 + i * 10 + 3);

        auto to_x_y_xyz_0 = pbuffer.data(idx_op_geom_101_fs + 1 * op_comps * 10 + i * 10 + 4);

        auto to_x_y_xzz_0 = pbuffer.data(idx_op_geom_101_fs + 1 * op_comps * 10 + i * 10 + 5);

        auto to_x_y_yyy_0 = pbuffer.data(idx_op_geom_101_fs + 1 * op_comps * 10 + i * 10 + 6);

        auto to_x_y_yyz_0 = pbuffer.data(idx_op_geom_101_fs + 1 * op_comps * 10 + i * 10 + 7);

        auto to_x_y_yzz_0 = pbuffer.data(idx_op_geom_101_fs + 1 * op_comps * 10 + i * 10 + 8);

        auto to_x_y_zzz_0 = pbuffer.data(idx_op_geom_101_fs + 1 * op_comps * 10 + i * 10 + 9);

#pragma omp simd aligned(to_x_y_xxx_0,     \
                             to_x_y_xxy_0, \
                             to_x_y_xxz_0, \
                             to_x_y_xyy_0, \
                             to_x_y_xyz_0, \
                             to_x_y_xzz_0, \
                             to_x_y_yyy_0, \
                             to_x_y_yyz_0, \
                             to_x_y_yzz_0, \
                             to_x_y_zzz_0, \
                             to_xx_y,      \
                             to_xxxx_y,    \
                             to_xxxy_y,    \
                             to_xxxz_y,    \
                             to_xxyy_y,    \
                             to_xxyz_y,    \
                             to_xxzz_y,    \
                             to_xy_y,      \
                             to_xyyy_y,    \
                             to_xyyz_y,    \
                             to_xyzz_y,    \
                             to_xz_y,      \
                             to_xzzz_y,    \
                             to_yy_y,      \
                             to_yz_y,      \
                             to_zz_y,      \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxx_0[k] = -6.0 * to_xx_y[k] * tke_0 + 4.0 * to_xxxx_y[k] * tbe_0 * tke_0;

            to_x_y_xxy_0[k] = -4.0 * to_xy_y[k] * tke_0 + 4.0 * to_xxxy_y[k] * tbe_0 * tke_0;

            to_x_y_xxz_0[k] = -4.0 * to_xz_y[k] * tke_0 + 4.0 * to_xxxz_y[k] * tbe_0 * tke_0;

            to_x_y_xyy_0[k] = -2.0 * to_yy_y[k] * tke_0 + 4.0 * to_xxyy_y[k] * tbe_0 * tke_0;

            to_x_y_xyz_0[k] = -2.0 * to_yz_y[k] * tke_0 + 4.0 * to_xxyz_y[k] * tbe_0 * tke_0;

            to_x_y_xzz_0[k] = -2.0 * to_zz_y[k] * tke_0 + 4.0 * to_xxzz_y[k] * tbe_0 * tke_0;

            to_x_y_yyy_0[k] = 4.0 * to_xyyy_y[k] * tbe_0 * tke_0;

            to_x_y_yyz_0[k] = 4.0 * to_xyyz_y[k] * tbe_0 * tke_0;

            to_x_y_yzz_0[k] = 4.0 * to_xyzz_y[k] * tbe_0 * tke_0;

            to_x_y_zzz_0[k] = 4.0 * to_xzzz_y[k] * tbe_0 * tke_0;
        }

        // Set up 20-30 components of targeted buffer : FS

        auto to_x_z_xxx_0 = pbuffer.data(idx_op_geom_101_fs + 2 * op_comps * 10 + i * 10 + 0);

        auto to_x_z_xxy_0 = pbuffer.data(idx_op_geom_101_fs + 2 * op_comps * 10 + i * 10 + 1);

        auto to_x_z_xxz_0 = pbuffer.data(idx_op_geom_101_fs + 2 * op_comps * 10 + i * 10 + 2);

        auto to_x_z_xyy_0 = pbuffer.data(idx_op_geom_101_fs + 2 * op_comps * 10 + i * 10 + 3);

        auto to_x_z_xyz_0 = pbuffer.data(idx_op_geom_101_fs + 2 * op_comps * 10 + i * 10 + 4);

        auto to_x_z_xzz_0 = pbuffer.data(idx_op_geom_101_fs + 2 * op_comps * 10 + i * 10 + 5);

        auto to_x_z_yyy_0 = pbuffer.data(idx_op_geom_101_fs + 2 * op_comps * 10 + i * 10 + 6);

        auto to_x_z_yyz_0 = pbuffer.data(idx_op_geom_101_fs + 2 * op_comps * 10 + i * 10 + 7);

        auto to_x_z_yzz_0 = pbuffer.data(idx_op_geom_101_fs + 2 * op_comps * 10 + i * 10 + 8);

        auto to_x_z_zzz_0 = pbuffer.data(idx_op_geom_101_fs + 2 * op_comps * 10 + i * 10 + 9);

#pragma omp simd aligned(to_x_z_xxx_0,     \
                             to_x_z_xxy_0, \
                             to_x_z_xxz_0, \
                             to_x_z_xyy_0, \
                             to_x_z_xyz_0, \
                             to_x_z_xzz_0, \
                             to_x_z_yyy_0, \
                             to_x_z_yyz_0, \
                             to_x_z_yzz_0, \
                             to_x_z_zzz_0, \
                             to_xx_z,      \
                             to_xxxx_z,    \
                             to_xxxy_z,    \
                             to_xxxz_z,    \
                             to_xxyy_z,    \
                             to_xxyz_z,    \
                             to_xxzz_z,    \
                             to_xy_z,      \
                             to_xyyy_z,    \
                             to_xyyz_z,    \
                             to_xyzz_z,    \
                             to_xz_z,      \
                             to_xzzz_z,    \
                             to_yy_z,      \
                             to_yz_z,      \
                             to_zz_z,      \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxx_0[k] = -6.0 * to_xx_z[k] * tke_0 + 4.0 * to_xxxx_z[k] * tbe_0 * tke_0;

            to_x_z_xxy_0[k] = -4.0 * to_xy_z[k] * tke_0 + 4.0 * to_xxxy_z[k] * tbe_0 * tke_0;

            to_x_z_xxz_0[k] = -4.0 * to_xz_z[k] * tke_0 + 4.0 * to_xxxz_z[k] * tbe_0 * tke_0;

            to_x_z_xyy_0[k] = -2.0 * to_yy_z[k] * tke_0 + 4.0 * to_xxyy_z[k] * tbe_0 * tke_0;

            to_x_z_xyz_0[k] = -2.0 * to_yz_z[k] * tke_0 + 4.0 * to_xxyz_z[k] * tbe_0 * tke_0;

            to_x_z_xzz_0[k] = -2.0 * to_zz_z[k] * tke_0 + 4.0 * to_xxzz_z[k] * tbe_0 * tke_0;

            to_x_z_yyy_0[k] = 4.0 * to_xyyy_z[k] * tbe_0 * tke_0;

            to_x_z_yyz_0[k] = 4.0 * to_xyyz_z[k] * tbe_0 * tke_0;

            to_x_z_yzz_0[k] = 4.0 * to_xyzz_z[k] * tbe_0 * tke_0;

            to_x_z_zzz_0[k] = 4.0 * to_xzzz_z[k] * tbe_0 * tke_0;
        }

        // Set up 30-40 components of targeted buffer : FS

        auto to_y_x_xxx_0 = pbuffer.data(idx_op_geom_101_fs + 3 * op_comps * 10 + i * 10 + 0);

        auto to_y_x_xxy_0 = pbuffer.data(idx_op_geom_101_fs + 3 * op_comps * 10 + i * 10 + 1);

        auto to_y_x_xxz_0 = pbuffer.data(idx_op_geom_101_fs + 3 * op_comps * 10 + i * 10 + 2);

        auto to_y_x_xyy_0 = pbuffer.data(idx_op_geom_101_fs + 3 * op_comps * 10 + i * 10 + 3);

        auto to_y_x_xyz_0 = pbuffer.data(idx_op_geom_101_fs + 3 * op_comps * 10 + i * 10 + 4);

        auto to_y_x_xzz_0 = pbuffer.data(idx_op_geom_101_fs + 3 * op_comps * 10 + i * 10 + 5);

        auto to_y_x_yyy_0 = pbuffer.data(idx_op_geom_101_fs + 3 * op_comps * 10 + i * 10 + 6);

        auto to_y_x_yyz_0 = pbuffer.data(idx_op_geom_101_fs + 3 * op_comps * 10 + i * 10 + 7);

        auto to_y_x_yzz_0 = pbuffer.data(idx_op_geom_101_fs + 3 * op_comps * 10 + i * 10 + 8);

        auto to_y_x_zzz_0 = pbuffer.data(idx_op_geom_101_fs + 3 * op_comps * 10 + i * 10 + 9);

#pragma omp simd aligned(to_xx_x,          \
                             to_xxxy_x,    \
                             to_xxyy_x,    \
                             to_xxyz_x,    \
                             to_xy_x,      \
                             to_xyyy_x,    \
                             to_xyyz_x,    \
                             to_xyzz_x,    \
                             to_xz_x,      \
                             to_y_x_xxx_0, \
                             to_y_x_xxy_0, \
                             to_y_x_xxz_0, \
                             to_y_x_xyy_0, \
                             to_y_x_xyz_0, \
                             to_y_x_xzz_0, \
                             to_y_x_yyy_0, \
                             to_y_x_yyz_0, \
                             to_y_x_yzz_0, \
                             to_y_x_zzz_0, \
                             to_yy_x,      \
                             to_yyyy_x,    \
                             to_yyyz_x,    \
                             to_yyzz_x,    \
                             to_yz_x,      \
                             to_yzzz_x,    \
                             to_zz_x,      \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxx_0[k] = 4.0 * to_xxxy_x[k] * tbe_0 * tke_0;

            to_y_x_xxy_0[k] = -2.0 * to_xx_x[k] * tke_0 + 4.0 * to_xxyy_x[k] * tbe_0 * tke_0;

            to_y_x_xxz_0[k] = 4.0 * to_xxyz_x[k] * tbe_0 * tke_0;

            to_y_x_xyy_0[k] = -4.0 * to_xy_x[k] * tke_0 + 4.0 * to_xyyy_x[k] * tbe_0 * tke_0;

            to_y_x_xyz_0[k] = -2.0 * to_xz_x[k] * tke_0 + 4.0 * to_xyyz_x[k] * tbe_0 * tke_0;

            to_y_x_xzz_0[k] = 4.0 * to_xyzz_x[k] * tbe_0 * tke_0;

            to_y_x_yyy_0[k] = -6.0 * to_yy_x[k] * tke_0 + 4.0 * to_yyyy_x[k] * tbe_0 * tke_0;

            to_y_x_yyz_0[k] = -4.0 * to_yz_x[k] * tke_0 + 4.0 * to_yyyz_x[k] * tbe_0 * tke_0;

            to_y_x_yzz_0[k] = -2.0 * to_zz_x[k] * tke_0 + 4.0 * to_yyzz_x[k] * tbe_0 * tke_0;

            to_y_x_zzz_0[k] = 4.0 * to_yzzz_x[k] * tbe_0 * tke_0;
        }

        // Set up 40-50 components of targeted buffer : FS

        auto to_y_y_xxx_0 = pbuffer.data(idx_op_geom_101_fs + 4 * op_comps * 10 + i * 10 + 0);

        auto to_y_y_xxy_0 = pbuffer.data(idx_op_geom_101_fs + 4 * op_comps * 10 + i * 10 + 1);

        auto to_y_y_xxz_0 = pbuffer.data(idx_op_geom_101_fs + 4 * op_comps * 10 + i * 10 + 2);

        auto to_y_y_xyy_0 = pbuffer.data(idx_op_geom_101_fs + 4 * op_comps * 10 + i * 10 + 3);

        auto to_y_y_xyz_0 = pbuffer.data(idx_op_geom_101_fs + 4 * op_comps * 10 + i * 10 + 4);

        auto to_y_y_xzz_0 = pbuffer.data(idx_op_geom_101_fs + 4 * op_comps * 10 + i * 10 + 5);

        auto to_y_y_yyy_0 = pbuffer.data(idx_op_geom_101_fs + 4 * op_comps * 10 + i * 10 + 6);

        auto to_y_y_yyz_0 = pbuffer.data(idx_op_geom_101_fs + 4 * op_comps * 10 + i * 10 + 7);

        auto to_y_y_yzz_0 = pbuffer.data(idx_op_geom_101_fs + 4 * op_comps * 10 + i * 10 + 8);

        auto to_y_y_zzz_0 = pbuffer.data(idx_op_geom_101_fs + 4 * op_comps * 10 + i * 10 + 9);

#pragma omp simd aligned(to_xx_y,          \
                             to_xxxy_y,    \
                             to_xxyy_y,    \
                             to_xxyz_y,    \
                             to_xy_y,      \
                             to_xyyy_y,    \
                             to_xyyz_y,    \
                             to_xyzz_y,    \
                             to_xz_y,      \
                             to_y_y_xxx_0, \
                             to_y_y_xxy_0, \
                             to_y_y_xxz_0, \
                             to_y_y_xyy_0, \
                             to_y_y_xyz_0, \
                             to_y_y_xzz_0, \
                             to_y_y_yyy_0, \
                             to_y_y_yyz_0, \
                             to_y_y_yzz_0, \
                             to_y_y_zzz_0, \
                             to_yy_y,      \
                             to_yyyy_y,    \
                             to_yyyz_y,    \
                             to_yyzz_y,    \
                             to_yz_y,      \
                             to_yzzz_y,    \
                             to_zz_y,      \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxx_0[k] = 4.0 * to_xxxy_y[k] * tbe_0 * tke_0;

            to_y_y_xxy_0[k] = -2.0 * to_xx_y[k] * tke_0 + 4.0 * to_xxyy_y[k] * tbe_0 * tke_0;

            to_y_y_xxz_0[k] = 4.0 * to_xxyz_y[k] * tbe_0 * tke_0;

            to_y_y_xyy_0[k] = -4.0 * to_xy_y[k] * tke_0 + 4.0 * to_xyyy_y[k] * tbe_0 * tke_0;

            to_y_y_xyz_0[k] = -2.0 * to_xz_y[k] * tke_0 + 4.0 * to_xyyz_y[k] * tbe_0 * tke_0;

            to_y_y_xzz_0[k] = 4.0 * to_xyzz_y[k] * tbe_0 * tke_0;

            to_y_y_yyy_0[k] = -6.0 * to_yy_y[k] * tke_0 + 4.0 * to_yyyy_y[k] * tbe_0 * tke_0;

            to_y_y_yyz_0[k] = -4.0 * to_yz_y[k] * tke_0 + 4.0 * to_yyyz_y[k] * tbe_0 * tke_0;

            to_y_y_yzz_0[k] = -2.0 * to_zz_y[k] * tke_0 + 4.0 * to_yyzz_y[k] * tbe_0 * tke_0;

            to_y_y_zzz_0[k] = 4.0 * to_yzzz_y[k] * tbe_0 * tke_0;
        }

        // Set up 50-60 components of targeted buffer : FS

        auto to_y_z_xxx_0 = pbuffer.data(idx_op_geom_101_fs + 5 * op_comps * 10 + i * 10 + 0);

        auto to_y_z_xxy_0 = pbuffer.data(idx_op_geom_101_fs + 5 * op_comps * 10 + i * 10 + 1);

        auto to_y_z_xxz_0 = pbuffer.data(idx_op_geom_101_fs + 5 * op_comps * 10 + i * 10 + 2);

        auto to_y_z_xyy_0 = pbuffer.data(idx_op_geom_101_fs + 5 * op_comps * 10 + i * 10 + 3);

        auto to_y_z_xyz_0 = pbuffer.data(idx_op_geom_101_fs + 5 * op_comps * 10 + i * 10 + 4);

        auto to_y_z_xzz_0 = pbuffer.data(idx_op_geom_101_fs + 5 * op_comps * 10 + i * 10 + 5);

        auto to_y_z_yyy_0 = pbuffer.data(idx_op_geom_101_fs + 5 * op_comps * 10 + i * 10 + 6);

        auto to_y_z_yyz_0 = pbuffer.data(idx_op_geom_101_fs + 5 * op_comps * 10 + i * 10 + 7);

        auto to_y_z_yzz_0 = pbuffer.data(idx_op_geom_101_fs + 5 * op_comps * 10 + i * 10 + 8);

        auto to_y_z_zzz_0 = pbuffer.data(idx_op_geom_101_fs + 5 * op_comps * 10 + i * 10 + 9);

#pragma omp simd aligned(to_xx_z,          \
                             to_xxxy_z,    \
                             to_xxyy_z,    \
                             to_xxyz_z,    \
                             to_xy_z,      \
                             to_xyyy_z,    \
                             to_xyyz_z,    \
                             to_xyzz_z,    \
                             to_xz_z,      \
                             to_y_z_xxx_0, \
                             to_y_z_xxy_0, \
                             to_y_z_xxz_0, \
                             to_y_z_xyy_0, \
                             to_y_z_xyz_0, \
                             to_y_z_xzz_0, \
                             to_y_z_yyy_0, \
                             to_y_z_yyz_0, \
                             to_y_z_yzz_0, \
                             to_y_z_zzz_0, \
                             to_yy_z,      \
                             to_yyyy_z,    \
                             to_yyyz_z,    \
                             to_yyzz_z,    \
                             to_yz_z,      \
                             to_yzzz_z,    \
                             to_zz_z,      \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxx_0[k] = 4.0 * to_xxxy_z[k] * tbe_0 * tke_0;

            to_y_z_xxy_0[k] = -2.0 * to_xx_z[k] * tke_0 + 4.0 * to_xxyy_z[k] * tbe_0 * tke_0;

            to_y_z_xxz_0[k] = 4.0 * to_xxyz_z[k] * tbe_0 * tke_0;

            to_y_z_xyy_0[k] = -4.0 * to_xy_z[k] * tke_0 + 4.0 * to_xyyy_z[k] * tbe_0 * tke_0;

            to_y_z_xyz_0[k] = -2.0 * to_xz_z[k] * tke_0 + 4.0 * to_xyyz_z[k] * tbe_0 * tke_0;

            to_y_z_xzz_0[k] = 4.0 * to_xyzz_z[k] * tbe_0 * tke_0;

            to_y_z_yyy_0[k] = -6.0 * to_yy_z[k] * tke_0 + 4.0 * to_yyyy_z[k] * tbe_0 * tke_0;

            to_y_z_yyz_0[k] = -4.0 * to_yz_z[k] * tke_0 + 4.0 * to_yyyz_z[k] * tbe_0 * tke_0;

            to_y_z_yzz_0[k] = -2.0 * to_zz_z[k] * tke_0 + 4.0 * to_yyzz_z[k] * tbe_0 * tke_0;

            to_y_z_zzz_0[k] = 4.0 * to_yzzz_z[k] * tbe_0 * tke_0;
        }

        // Set up 60-70 components of targeted buffer : FS

        auto to_z_x_xxx_0 = pbuffer.data(idx_op_geom_101_fs + 6 * op_comps * 10 + i * 10 + 0);

        auto to_z_x_xxy_0 = pbuffer.data(idx_op_geom_101_fs + 6 * op_comps * 10 + i * 10 + 1);

        auto to_z_x_xxz_0 = pbuffer.data(idx_op_geom_101_fs + 6 * op_comps * 10 + i * 10 + 2);

        auto to_z_x_xyy_0 = pbuffer.data(idx_op_geom_101_fs + 6 * op_comps * 10 + i * 10 + 3);

        auto to_z_x_xyz_0 = pbuffer.data(idx_op_geom_101_fs + 6 * op_comps * 10 + i * 10 + 4);

        auto to_z_x_xzz_0 = pbuffer.data(idx_op_geom_101_fs + 6 * op_comps * 10 + i * 10 + 5);

        auto to_z_x_yyy_0 = pbuffer.data(idx_op_geom_101_fs + 6 * op_comps * 10 + i * 10 + 6);

        auto to_z_x_yyz_0 = pbuffer.data(idx_op_geom_101_fs + 6 * op_comps * 10 + i * 10 + 7);

        auto to_z_x_yzz_0 = pbuffer.data(idx_op_geom_101_fs + 6 * op_comps * 10 + i * 10 + 8);

        auto to_z_x_zzz_0 = pbuffer.data(idx_op_geom_101_fs + 6 * op_comps * 10 + i * 10 + 9);

#pragma omp simd aligned(to_xx_x,          \
                             to_xxxz_x,    \
                             to_xxyz_x,    \
                             to_xxzz_x,    \
                             to_xy_x,      \
                             to_xyyz_x,    \
                             to_xyzz_x,    \
                             to_xz_x,      \
                             to_xzzz_x,    \
                             to_yy_x,      \
                             to_yyyz_x,    \
                             to_yyzz_x,    \
                             to_yz_x,      \
                             to_yzzz_x,    \
                             to_z_x_xxx_0, \
                             to_z_x_xxy_0, \
                             to_z_x_xxz_0, \
                             to_z_x_xyy_0, \
                             to_z_x_xyz_0, \
                             to_z_x_xzz_0, \
                             to_z_x_yyy_0, \
                             to_z_x_yyz_0, \
                             to_z_x_yzz_0, \
                             to_z_x_zzz_0, \
                             to_zz_x,      \
                             to_zzzz_x,    \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxx_0[k] = 4.0 * to_xxxz_x[k] * tbe_0 * tke_0;

            to_z_x_xxy_0[k] = 4.0 * to_xxyz_x[k] * tbe_0 * tke_0;

            to_z_x_xxz_0[k] = -2.0 * to_xx_x[k] * tke_0 + 4.0 * to_xxzz_x[k] * tbe_0 * tke_0;

            to_z_x_xyy_0[k] = 4.0 * to_xyyz_x[k] * tbe_0 * tke_0;

            to_z_x_xyz_0[k] = -2.0 * to_xy_x[k] * tke_0 + 4.0 * to_xyzz_x[k] * tbe_0 * tke_0;

            to_z_x_xzz_0[k] = -4.0 * to_xz_x[k] * tke_0 + 4.0 * to_xzzz_x[k] * tbe_0 * tke_0;

            to_z_x_yyy_0[k] = 4.0 * to_yyyz_x[k] * tbe_0 * tke_0;

            to_z_x_yyz_0[k] = -2.0 * to_yy_x[k] * tke_0 + 4.0 * to_yyzz_x[k] * tbe_0 * tke_0;

            to_z_x_yzz_0[k] = -4.0 * to_yz_x[k] * tke_0 + 4.0 * to_yzzz_x[k] * tbe_0 * tke_0;

            to_z_x_zzz_0[k] = -6.0 * to_zz_x[k] * tke_0 + 4.0 * to_zzzz_x[k] * tbe_0 * tke_0;
        }

        // Set up 70-80 components of targeted buffer : FS

        auto to_z_y_xxx_0 = pbuffer.data(idx_op_geom_101_fs + 7 * op_comps * 10 + i * 10 + 0);

        auto to_z_y_xxy_0 = pbuffer.data(idx_op_geom_101_fs + 7 * op_comps * 10 + i * 10 + 1);

        auto to_z_y_xxz_0 = pbuffer.data(idx_op_geom_101_fs + 7 * op_comps * 10 + i * 10 + 2);

        auto to_z_y_xyy_0 = pbuffer.data(idx_op_geom_101_fs + 7 * op_comps * 10 + i * 10 + 3);

        auto to_z_y_xyz_0 = pbuffer.data(idx_op_geom_101_fs + 7 * op_comps * 10 + i * 10 + 4);

        auto to_z_y_xzz_0 = pbuffer.data(idx_op_geom_101_fs + 7 * op_comps * 10 + i * 10 + 5);

        auto to_z_y_yyy_0 = pbuffer.data(idx_op_geom_101_fs + 7 * op_comps * 10 + i * 10 + 6);

        auto to_z_y_yyz_0 = pbuffer.data(idx_op_geom_101_fs + 7 * op_comps * 10 + i * 10 + 7);

        auto to_z_y_yzz_0 = pbuffer.data(idx_op_geom_101_fs + 7 * op_comps * 10 + i * 10 + 8);

        auto to_z_y_zzz_0 = pbuffer.data(idx_op_geom_101_fs + 7 * op_comps * 10 + i * 10 + 9);

#pragma omp simd aligned(to_xx_y,          \
                             to_xxxz_y,    \
                             to_xxyz_y,    \
                             to_xxzz_y,    \
                             to_xy_y,      \
                             to_xyyz_y,    \
                             to_xyzz_y,    \
                             to_xz_y,      \
                             to_xzzz_y,    \
                             to_yy_y,      \
                             to_yyyz_y,    \
                             to_yyzz_y,    \
                             to_yz_y,      \
                             to_yzzz_y,    \
                             to_z_y_xxx_0, \
                             to_z_y_xxy_0, \
                             to_z_y_xxz_0, \
                             to_z_y_xyy_0, \
                             to_z_y_xyz_0, \
                             to_z_y_xzz_0, \
                             to_z_y_yyy_0, \
                             to_z_y_yyz_0, \
                             to_z_y_yzz_0, \
                             to_z_y_zzz_0, \
                             to_zz_y,      \
                             to_zzzz_y,    \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxx_0[k] = 4.0 * to_xxxz_y[k] * tbe_0 * tke_0;

            to_z_y_xxy_0[k] = 4.0 * to_xxyz_y[k] * tbe_0 * tke_0;

            to_z_y_xxz_0[k] = -2.0 * to_xx_y[k] * tke_0 + 4.0 * to_xxzz_y[k] * tbe_0 * tke_0;

            to_z_y_xyy_0[k] = 4.0 * to_xyyz_y[k] * tbe_0 * tke_0;

            to_z_y_xyz_0[k] = -2.0 * to_xy_y[k] * tke_0 + 4.0 * to_xyzz_y[k] * tbe_0 * tke_0;

            to_z_y_xzz_0[k] = -4.0 * to_xz_y[k] * tke_0 + 4.0 * to_xzzz_y[k] * tbe_0 * tke_0;

            to_z_y_yyy_0[k] = 4.0 * to_yyyz_y[k] * tbe_0 * tke_0;

            to_z_y_yyz_0[k] = -2.0 * to_yy_y[k] * tke_0 + 4.0 * to_yyzz_y[k] * tbe_0 * tke_0;

            to_z_y_yzz_0[k] = -4.0 * to_yz_y[k] * tke_0 + 4.0 * to_yzzz_y[k] * tbe_0 * tke_0;

            to_z_y_zzz_0[k] = -6.0 * to_zz_y[k] * tke_0 + 4.0 * to_zzzz_y[k] * tbe_0 * tke_0;
        }

        // Set up 80-90 components of targeted buffer : FS

        auto to_z_z_xxx_0 = pbuffer.data(idx_op_geom_101_fs + 8 * op_comps * 10 + i * 10 + 0);

        auto to_z_z_xxy_0 = pbuffer.data(idx_op_geom_101_fs + 8 * op_comps * 10 + i * 10 + 1);

        auto to_z_z_xxz_0 = pbuffer.data(idx_op_geom_101_fs + 8 * op_comps * 10 + i * 10 + 2);

        auto to_z_z_xyy_0 = pbuffer.data(idx_op_geom_101_fs + 8 * op_comps * 10 + i * 10 + 3);

        auto to_z_z_xyz_0 = pbuffer.data(idx_op_geom_101_fs + 8 * op_comps * 10 + i * 10 + 4);

        auto to_z_z_xzz_0 = pbuffer.data(idx_op_geom_101_fs + 8 * op_comps * 10 + i * 10 + 5);

        auto to_z_z_yyy_0 = pbuffer.data(idx_op_geom_101_fs + 8 * op_comps * 10 + i * 10 + 6);

        auto to_z_z_yyz_0 = pbuffer.data(idx_op_geom_101_fs + 8 * op_comps * 10 + i * 10 + 7);

        auto to_z_z_yzz_0 = pbuffer.data(idx_op_geom_101_fs + 8 * op_comps * 10 + i * 10 + 8);

        auto to_z_z_zzz_0 = pbuffer.data(idx_op_geom_101_fs + 8 * op_comps * 10 + i * 10 + 9);

#pragma omp simd aligned(to_xx_z,          \
                             to_xxxz_z,    \
                             to_xxyz_z,    \
                             to_xxzz_z,    \
                             to_xy_z,      \
                             to_xyyz_z,    \
                             to_xyzz_z,    \
                             to_xz_z,      \
                             to_xzzz_z,    \
                             to_yy_z,      \
                             to_yyyz_z,    \
                             to_yyzz_z,    \
                             to_yz_z,      \
                             to_yzzz_z,    \
                             to_z_z_xxx_0, \
                             to_z_z_xxy_0, \
                             to_z_z_xxz_0, \
                             to_z_z_xyy_0, \
                             to_z_z_xyz_0, \
                             to_z_z_xzz_0, \
                             to_z_z_yyy_0, \
                             to_z_z_yyz_0, \
                             to_z_z_yzz_0, \
                             to_z_z_zzz_0, \
                             to_zz_z,      \
                             to_zzzz_z,    \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxx_0[k] = 4.0 * to_xxxz_z[k] * tbe_0 * tke_0;

            to_z_z_xxy_0[k] = 4.0 * to_xxyz_z[k] * tbe_0 * tke_0;

            to_z_z_xxz_0[k] = -2.0 * to_xx_z[k] * tke_0 + 4.0 * to_xxzz_z[k] * tbe_0 * tke_0;

            to_z_z_xyy_0[k] = 4.0 * to_xyyz_z[k] * tbe_0 * tke_0;

            to_z_z_xyz_0[k] = -2.0 * to_xy_z[k] * tke_0 + 4.0 * to_xyzz_z[k] * tbe_0 * tke_0;

            to_z_z_xzz_0[k] = -4.0 * to_xz_z[k] * tke_0 + 4.0 * to_xzzz_z[k] * tbe_0 * tke_0;

            to_z_z_yyy_0[k] = 4.0 * to_yyyz_z[k] * tbe_0 * tke_0;

            to_z_z_yyz_0[k] = -2.0 * to_yy_z[k] * tke_0 + 4.0 * to_yyzz_z[k] * tbe_0 * tke_0;

            to_z_z_yzz_0[k] = -4.0 * to_yz_z[k] * tke_0 + 4.0 * to_yzzz_z[k] * tbe_0 * tke_0;

            to_z_z_zzz_0[k] = -6.0 * to_zz_z[k] * tke_0 + 4.0 * to_zzzz_z[k] * tbe_0 * tke_0;
        }
    }
}

}  // namespace t2cgeom
