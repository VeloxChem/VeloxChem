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

#include "GeometricalDerivatives1X1ForDS.hpp"

namespace t2cgeom {  // t2cgeom namespace

auto
comp_prim_op_geom_11_ds(CSimdArray<double>&       pbuffer,
                        const size_t              idx_op_geom_101_ds,
                        const size_t              idx_op_pp,
                        const size_t              idx_op_fp,
                        const size_t              op_comps,
                        const CSimdArray<double>& factors,
                        const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : PP

        auto to_x_x = pbuffer.data(idx_op_pp + i * 9 + 0);

        auto to_x_y = pbuffer.data(idx_op_pp + i * 9 + 1);

        auto to_x_z = pbuffer.data(idx_op_pp + i * 9 + 2);

        auto to_y_x = pbuffer.data(idx_op_pp + i * 9 + 3);

        auto to_y_y = pbuffer.data(idx_op_pp + i * 9 + 4);

        auto to_y_z = pbuffer.data(idx_op_pp + i * 9 + 5);

        auto to_z_x = pbuffer.data(idx_op_pp + i * 9 + 6);

        auto to_z_y = pbuffer.data(idx_op_pp + i * 9 + 7);

        auto to_z_z = pbuffer.data(idx_op_pp + i * 9 + 8);

        // Set up components of auxiliary buffer : FP

        auto to_xxx_x = pbuffer.data(idx_op_fp + i * 30 + 0);

        auto to_xxx_y = pbuffer.data(idx_op_fp + i * 30 + 1);

        auto to_xxx_z = pbuffer.data(idx_op_fp + i * 30 + 2);

        auto to_xxy_x = pbuffer.data(idx_op_fp + i * 30 + 3);

        auto to_xxy_y = pbuffer.data(idx_op_fp + i * 30 + 4);

        auto to_xxy_z = pbuffer.data(idx_op_fp + i * 30 + 5);

        auto to_xxz_x = pbuffer.data(idx_op_fp + i * 30 + 6);

        auto to_xxz_y = pbuffer.data(idx_op_fp + i * 30 + 7);

        auto to_xxz_z = pbuffer.data(idx_op_fp + i * 30 + 8);

        auto to_xyy_x = pbuffer.data(idx_op_fp + i * 30 + 9);

        auto to_xyy_y = pbuffer.data(idx_op_fp + i * 30 + 10);

        auto to_xyy_z = pbuffer.data(idx_op_fp + i * 30 + 11);

        auto to_xyz_x = pbuffer.data(idx_op_fp + i * 30 + 12);

        auto to_xyz_y = pbuffer.data(idx_op_fp + i * 30 + 13);

        auto to_xyz_z = pbuffer.data(idx_op_fp + i * 30 + 14);

        auto to_xzz_x = pbuffer.data(idx_op_fp + i * 30 + 15);

        auto to_xzz_y = pbuffer.data(idx_op_fp + i * 30 + 16);

        auto to_xzz_z = pbuffer.data(idx_op_fp + i * 30 + 17);

        auto to_yyy_x = pbuffer.data(idx_op_fp + i * 30 + 18);

        auto to_yyy_y = pbuffer.data(idx_op_fp + i * 30 + 19);

        auto to_yyy_z = pbuffer.data(idx_op_fp + i * 30 + 20);

        auto to_yyz_x = pbuffer.data(idx_op_fp + i * 30 + 21);

        auto to_yyz_y = pbuffer.data(idx_op_fp + i * 30 + 22);

        auto to_yyz_z = pbuffer.data(idx_op_fp + i * 30 + 23);

        auto to_yzz_x = pbuffer.data(idx_op_fp + i * 30 + 24);

        auto to_yzz_y = pbuffer.data(idx_op_fp + i * 30 + 25);

        auto to_yzz_z = pbuffer.data(idx_op_fp + i * 30 + 26);

        auto to_zzz_x = pbuffer.data(idx_op_fp + i * 30 + 27);

        auto to_zzz_y = pbuffer.data(idx_op_fp + i * 30 + 28);

        auto to_zzz_z = pbuffer.data(idx_op_fp + i * 30 + 29);

        // Set up 0-6 components of targeted buffer : DS

        auto to_x_x_xx_0 = pbuffer.data(idx_op_geom_101_ds + 0 * op_comps * 6 + i * 6 + 0);

        auto to_x_x_xy_0 = pbuffer.data(idx_op_geom_101_ds + 0 * op_comps * 6 + i * 6 + 1);

        auto to_x_x_xz_0 = pbuffer.data(idx_op_geom_101_ds + 0 * op_comps * 6 + i * 6 + 2);

        auto to_x_x_yy_0 = pbuffer.data(idx_op_geom_101_ds + 0 * op_comps * 6 + i * 6 + 3);

        auto to_x_x_yz_0 = pbuffer.data(idx_op_geom_101_ds + 0 * op_comps * 6 + i * 6 + 4);

        auto to_x_x_zz_0 = pbuffer.data(idx_op_geom_101_ds + 0 * op_comps * 6 + i * 6 + 5);

#pragma omp simd aligned(to_x_x,          \
                             to_x_x_xx_0, \
                             to_x_x_xy_0, \
                             to_x_x_xz_0, \
                             to_x_x_yy_0, \
                             to_x_x_yz_0, \
                             to_x_x_zz_0, \
                             to_xxx_x,    \
                             to_xxy_x,    \
                             to_xxz_x,    \
                             to_xyy_x,    \
                             to_xyz_x,    \
                             to_xzz_x,    \
                             to_y_x,      \
                             to_z_x,      \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xx_0[k] = -4.0 * to_x_x[k] * tke_0 + 4.0 * to_xxx_x[k] * tbe_0 * tke_0;

            to_x_x_xy_0[k] = -2.0 * to_y_x[k] * tke_0 + 4.0 * to_xxy_x[k] * tbe_0 * tke_0;

            to_x_x_xz_0[k] = -2.0 * to_z_x[k] * tke_0 + 4.0 * to_xxz_x[k] * tbe_0 * tke_0;

            to_x_x_yy_0[k] = 4.0 * to_xyy_x[k] * tbe_0 * tke_0;

            to_x_x_yz_0[k] = 4.0 * to_xyz_x[k] * tbe_0 * tke_0;

            to_x_x_zz_0[k] = 4.0 * to_xzz_x[k] * tbe_0 * tke_0;
        }

        // Set up 6-12 components of targeted buffer : DS

        auto to_x_y_xx_0 = pbuffer.data(idx_op_geom_101_ds + 1 * op_comps * 6 + i * 6 + 0);

        auto to_x_y_xy_0 = pbuffer.data(idx_op_geom_101_ds + 1 * op_comps * 6 + i * 6 + 1);

        auto to_x_y_xz_0 = pbuffer.data(idx_op_geom_101_ds + 1 * op_comps * 6 + i * 6 + 2);

        auto to_x_y_yy_0 = pbuffer.data(idx_op_geom_101_ds + 1 * op_comps * 6 + i * 6 + 3);

        auto to_x_y_yz_0 = pbuffer.data(idx_op_geom_101_ds + 1 * op_comps * 6 + i * 6 + 4);

        auto to_x_y_zz_0 = pbuffer.data(idx_op_geom_101_ds + 1 * op_comps * 6 + i * 6 + 5);

#pragma omp simd aligned(to_x_y,          \
                             to_x_y_xx_0, \
                             to_x_y_xy_0, \
                             to_x_y_xz_0, \
                             to_x_y_yy_0, \
                             to_x_y_yz_0, \
                             to_x_y_zz_0, \
                             to_xxx_y,    \
                             to_xxy_y,    \
                             to_xxz_y,    \
                             to_xyy_y,    \
                             to_xyz_y,    \
                             to_xzz_y,    \
                             to_y_y,      \
                             to_z_y,      \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xx_0[k] = -4.0 * to_x_y[k] * tke_0 + 4.0 * to_xxx_y[k] * tbe_0 * tke_0;

            to_x_y_xy_0[k] = -2.0 * to_y_y[k] * tke_0 + 4.0 * to_xxy_y[k] * tbe_0 * tke_0;

            to_x_y_xz_0[k] = -2.0 * to_z_y[k] * tke_0 + 4.0 * to_xxz_y[k] * tbe_0 * tke_0;

            to_x_y_yy_0[k] = 4.0 * to_xyy_y[k] * tbe_0 * tke_0;

            to_x_y_yz_0[k] = 4.0 * to_xyz_y[k] * tbe_0 * tke_0;

            to_x_y_zz_0[k] = 4.0 * to_xzz_y[k] * tbe_0 * tke_0;
        }

        // Set up 12-18 components of targeted buffer : DS

        auto to_x_z_xx_0 = pbuffer.data(idx_op_geom_101_ds + 2 * op_comps * 6 + i * 6 + 0);

        auto to_x_z_xy_0 = pbuffer.data(idx_op_geom_101_ds + 2 * op_comps * 6 + i * 6 + 1);

        auto to_x_z_xz_0 = pbuffer.data(idx_op_geom_101_ds + 2 * op_comps * 6 + i * 6 + 2);

        auto to_x_z_yy_0 = pbuffer.data(idx_op_geom_101_ds + 2 * op_comps * 6 + i * 6 + 3);

        auto to_x_z_yz_0 = pbuffer.data(idx_op_geom_101_ds + 2 * op_comps * 6 + i * 6 + 4);

        auto to_x_z_zz_0 = pbuffer.data(idx_op_geom_101_ds + 2 * op_comps * 6 + i * 6 + 5);

#pragma omp simd aligned(to_x_z,          \
                             to_x_z_xx_0, \
                             to_x_z_xy_0, \
                             to_x_z_xz_0, \
                             to_x_z_yy_0, \
                             to_x_z_yz_0, \
                             to_x_z_zz_0, \
                             to_xxx_z,    \
                             to_xxy_z,    \
                             to_xxz_z,    \
                             to_xyy_z,    \
                             to_xyz_z,    \
                             to_xzz_z,    \
                             to_y_z,      \
                             to_z_z,      \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xx_0[k] = -4.0 * to_x_z[k] * tke_0 + 4.0 * to_xxx_z[k] * tbe_0 * tke_0;

            to_x_z_xy_0[k] = -2.0 * to_y_z[k] * tke_0 + 4.0 * to_xxy_z[k] * tbe_0 * tke_0;

            to_x_z_xz_0[k] = -2.0 * to_z_z[k] * tke_0 + 4.0 * to_xxz_z[k] * tbe_0 * tke_0;

            to_x_z_yy_0[k] = 4.0 * to_xyy_z[k] * tbe_0 * tke_0;

            to_x_z_yz_0[k] = 4.0 * to_xyz_z[k] * tbe_0 * tke_0;

            to_x_z_zz_0[k] = 4.0 * to_xzz_z[k] * tbe_0 * tke_0;
        }

        // Set up 18-24 components of targeted buffer : DS

        auto to_y_x_xx_0 = pbuffer.data(idx_op_geom_101_ds + 3 * op_comps * 6 + i * 6 + 0);

        auto to_y_x_xy_0 = pbuffer.data(idx_op_geom_101_ds + 3 * op_comps * 6 + i * 6 + 1);

        auto to_y_x_xz_0 = pbuffer.data(idx_op_geom_101_ds + 3 * op_comps * 6 + i * 6 + 2);

        auto to_y_x_yy_0 = pbuffer.data(idx_op_geom_101_ds + 3 * op_comps * 6 + i * 6 + 3);

        auto to_y_x_yz_0 = pbuffer.data(idx_op_geom_101_ds + 3 * op_comps * 6 + i * 6 + 4);

        auto to_y_x_zz_0 = pbuffer.data(idx_op_geom_101_ds + 3 * op_comps * 6 + i * 6 + 5);

#pragma omp simd aligned(to_x_x,          \
                             to_xxy_x,    \
                             to_xyy_x,    \
                             to_xyz_x,    \
                             to_y_x,      \
                             to_y_x_xx_0, \
                             to_y_x_xy_0, \
                             to_y_x_xz_0, \
                             to_y_x_yy_0, \
                             to_y_x_yz_0, \
                             to_y_x_zz_0, \
                             to_yyy_x,    \
                             to_yyz_x,    \
                             to_yzz_x,    \
                             to_z_x,      \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xx_0[k] = 4.0 * to_xxy_x[k] * tbe_0 * tke_0;

            to_y_x_xy_0[k] = -2.0 * to_x_x[k] * tke_0 + 4.0 * to_xyy_x[k] * tbe_0 * tke_0;

            to_y_x_xz_0[k] = 4.0 * to_xyz_x[k] * tbe_0 * tke_0;

            to_y_x_yy_0[k] = -4.0 * to_y_x[k] * tke_0 + 4.0 * to_yyy_x[k] * tbe_0 * tke_0;

            to_y_x_yz_0[k] = -2.0 * to_z_x[k] * tke_0 + 4.0 * to_yyz_x[k] * tbe_0 * tke_0;

            to_y_x_zz_0[k] = 4.0 * to_yzz_x[k] * tbe_0 * tke_0;
        }

        // Set up 24-30 components of targeted buffer : DS

        auto to_y_y_xx_0 = pbuffer.data(idx_op_geom_101_ds + 4 * op_comps * 6 + i * 6 + 0);

        auto to_y_y_xy_0 = pbuffer.data(idx_op_geom_101_ds + 4 * op_comps * 6 + i * 6 + 1);

        auto to_y_y_xz_0 = pbuffer.data(idx_op_geom_101_ds + 4 * op_comps * 6 + i * 6 + 2);

        auto to_y_y_yy_0 = pbuffer.data(idx_op_geom_101_ds + 4 * op_comps * 6 + i * 6 + 3);

        auto to_y_y_yz_0 = pbuffer.data(idx_op_geom_101_ds + 4 * op_comps * 6 + i * 6 + 4);

        auto to_y_y_zz_0 = pbuffer.data(idx_op_geom_101_ds + 4 * op_comps * 6 + i * 6 + 5);

#pragma omp simd aligned(to_x_y,          \
                             to_xxy_y,    \
                             to_xyy_y,    \
                             to_xyz_y,    \
                             to_y_y,      \
                             to_y_y_xx_0, \
                             to_y_y_xy_0, \
                             to_y_y_xz_0, \
                             to_y_y_yy_0, \
                             to_y_y_yz_0, \
                             to_y_y_zz_0, \
                             to_yyy_y,    \
                             to_yyz_y,    \
                             to_yzz_y,    \
                             to_z_y,      \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xx_0[k] = 4.0 * to_xxy_y[k] * tbe_0 * tke_0;

            to_y_y_xy_0[k] = -2.0 * to_x_y[k] * tke_0 + 4.0 * to_xyy_y[k] * tbe_0 * tke_0;

            to_y_y_xz_0[k] = 4.0 * to_xyz_y[k] * tbe_0 * tke_0;

            to_y_y_yy_0[k] = -4.0 * to_y_y[k] * tke_0 + 4.0 * to_yyy_y[k] * tbe_0 * tke_0;

            to_y_y_yz_0[k] = -2.0 * to_z_y[k] * tke_0 + 4.0 * to_yyz_y[k] * tbe_0 * tke_0;

            to_y_y_zz_0[k] = 4.0 * to_yzz_y[k] * tbe_0 * tke_0;
        }

        // Set up 30-36 components of targeted buffer : DS

        auto to_y_z_xx_0 = pbuffer.data(idx_op_geom_101_ds + 5 * op_comps * 6 + i * 6 + 0);

        auto to_y_z_xy_0 = pbuffer.data(idx_op_geom_101_ds + 5 * op_comps * 6 + i * 6 + 1);

        auto to_y_z_xz_0 = pbuffer.data(idx_op_geom_101_ds + 5 * op_comps * 6 + i * 6 + 2);

        auto to_y_z_yy_0 = pbuffer.data(idx_op_geom_101_ds + 5 * op_comps * 6 + i * 6 + 3);

        auto to_y_z_yz_0 = pbuffer.data(idx_op_geom_101_ds + 5 * op_comps * 6 + i * 6 + 4);

        auto to_y_z_zz_0 = pbuffer.data(idx_op_geom_101_ds + 5 * op_comps * 6 + i * 6 + 5);

#pragma omp simd aligned(to_x_z,          \
                             to_xxy_z,    \
                             to_xyy_z,    \
                             to_xyz_z,    \
                             to_y_z,      \
                             to_y_z_xx_0, \
                             to_y_z_xy_0, \
                             to_y_z_xz_0, \
                             to_y_z_yy_0, \
                             to_y_z_yz_0, \
                             to_y_z_zz_0, \
                             to_yyy_z,    \
                             to_yyz_z,    \
                             to_yzz_z,    \
                             to_z_z,      \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xx_0[k] = 4.0 * to_xxy_z[k] * tbe_0 * tke_0;

            to_y_z_xy_0[k] = -2.0 * to_x_z[k] * tke_0 + 4.0 * to_xyy_z[k] * tbe_0 * tke_0;

            to_y_z_xz_0[k] = 4.0 * to_xyz_z[k] * tbe_0 * tke_0;

            to_y_z_yy_0[k] = -4.0 * to_y_z[k] * tke_0 + 4.0 * to_yyy_z[k] * tbe_0 * tke_0;

            to_y_z_yz_0[k] = -2.0 * to_z_z[k] * tke_0 + 4.0 * to_yyz_z[k] * tbe_0 * tke_0;

            to_y_z_zz_0[k] = 4.0 * to_yzz_z[k] * tbe_0 * tke_0;
        }

        // Set up 36-42 components of targeted buffer : DS

        auto to_z_x_xx_0 = pbuffer.data(idx_op_geom_101_ds + 6 * op_comps * 6 + i * 6 + 0);

        auto to_z_x_xy_0 = pbuffer.data(idx_op_geom_101_ds + 6 * op_comps * 6 + i * 6 + 1);

        auto to_z_x_xz_0 = pbuffer.data(idx_op_geom_101_ds + 6 * op_comps * 6 + i * 6 + 2);

        auto to_z_x_yy_0 = pbuffer.data(idx_op_geom_101_ds + 6 * op_comps * 6 + i * 6 + 3);

        auto to_z_x_yz_0 = pbuffer.data(idx_op_geom_101_ds + 6 * op_comps * 6 + i * 6 + 4);

        auto to_z_x_zz_0 = pbuffer.data(idx_op_geom_101_ds + 6 * op_comps * 6 + i * 6 + 5);

#pragma omp simd aligned(to_x_x,          \
                             to_xxz_x,    \
                             to_xyz_x,    \
                             to_xzz_x,    \
                             to_y_x,      \
                             to_yyz_x,    \
                             to_yzz_x,    \
                             to_z_x,      \
                             to_z_x_xx_0, \
                             to_z_x_xy_0, \
                             to_z_x_xz_0, \
                             to_z_x_yy_0, \
                             to_z_x_yz_0, \
                             to_z_x_zz_0, \
                             to_zzz_x,    \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xx_0[k] = 4.0 * to_xxz_x[k] * tbe_0 * tke_0;

            to_z_x_xy_0[k] = 4.0 * to_xyz_x[k] * tbe_0 * tke_0;

            to_z_x_xz_0[k] = -2.0 * to_x_x[k] * tke_0 + 4.0 * to_xzz_x[k] * tbe_0 * tke_0;

            to_z_x_yy_0[k] = 4.0 * to_yyz_x[k] * tbe_0 * tke_0;

            to_z_x_yz_0[k] = -2.0 * to_y_x[k] * tke_0 + 4.0 * to_yzz_x[k] * tbe_0 * tke_0;

            to_z_x_zz_0[k] = -4.0 * to_z_x[k] * tke_0 + 4.0 * to_zzz_x[k] * tbe_0 * tke_0;
        }

        // Set up 42-48 components of targeted buffer : DS

        auto to_z_y_xx_0 = pbuffer.data(idx_op_geom_101_ds + 7 * op_comps * 6 + i * 6 + 0);

        auto to_z_y_xy_0 = pbuffer.data(idx_op_geom_101_ds + 7 * op_comps * 6 + i * 6 + 1);

        auto to_z_y_xz_0 = pbuffer.data(idx_op_geom_101_ds + 7 * op_comps * 6 + i * 6 + 2);

        auto to_z_y_yy_0 = pbuffer.data(idx_op_geom_101_ds + 7 * op_comps * 6 + i * 6 + 3);

        auto to_z_y_yz_0 = pbuffer.data(idx_op_geom_101_ds + 7 * op_comps * 6 + i * 6 + 4);

        auto to_z_y_zz_0 = pbuffer.data(idx_op_geom_101_ds + 7 * op_comps * 6 + i * 6 + 5);

#pragma omp simd aligned(to_x_y,          \
                             to_xxz_y,    \
                             to_xyz_y,    \
                             to_xzz_y,    \
                             to_y_y,      \
                             to_yyz_y,    \
                             to_yzz_y,    \
                             to_z_y,      \
                             to_z_y_xx_0, \
                             to_z_y_xy_0, \
                             to_z_y_xz_0, \
                             to_z_y_yy_0, \
                             to_z_y_yz_0, \
                             to_z_y_zz_0, \
                             to_zzz_y,    \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xx_0[k] = 4.0 * to_xxz_y[k] * tbe_0 * tke_0;

            to_z_y_xy_0[k] = 4.0 * to_xyz_y[k] * tbe_0 * tke_0;

            to_z_y_xz_0[k] = -2.0 * to_x_y[k] * tke_0 + 4.0 * to_xzz_y[k] * tbe_0 * tke_0;

            to_z_y_yy_0[k] = 4.0 * to_yyz_y[k] * tbe_0 * tke_0;

            to_z_y_yz_0[k] = -2.0 * to_y_y[k] * tke_0 + 4.0 * to_yzz_y[k] * tbe_0 * tke_0;

            to_z_y_zz_0[k] = -4.0 * to_z_y[k] * tke_0 + 4.0 * to_zzz_y[k] * tbe_0 * tke_0;
        }

        // Set up 48-54 components of targeted buffer : DS

        auto to_z_z_xx_0 = pbuffer.data(idx_op_geom_101_ds + 8 * op_comps * 6 + i * 6 + 0);

        auto to_z_z_xy_0 = pbuffer.data(idx_op_geom_101_ds + 8 * op_comps * 6 + i * 6 + 1);

        auto to_z_z_xz_0 = pbuffer.data(idx_op_geom_101_ds + 8 * op_comps * 6 + i * 6 + 2);

        auto to_z_z_yy_0 = pbuffer.data(idx_op_geom_101_ds + 8 * op_comps * 6 + i * 6 + 3);

        auto to_z_z_yz_0 = pbuffer.data(idx_op_geom_101_ds + 8 * op_comps * 6 + i * 6 + 4);

        auto to_z_z_zz_0 = pbuffer.data(idx_op_geom_101_ds + 8 * op_comps * 6 + i * 6 + 5);

#pragma omp simd aligned(to_x_z,          \
                             to_xxz_z,    \
                             to_xyz_z,    \
                             to_xzz_z,    \
                             to_y_z,      \
                             to_yyz_z,    \
                             to_yzz_z,    \
                             to_z_z,      \
                             to_z_z_xx_0, \
                             to_z_z_xy_0, \
                             to_z_z_xz_0, \
                             to_z_z_yy_0, \
                             to_z_z_yz_0, \
                             to_z_z_zz_0, \
                             to_zzz_z,    \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xx_0[k] = 4.0 * to_xxz_z[k] * tbe_0 * tke_0;

            to_z_z_xy_0[k] = 4.0 * to_xyz_z[k] * tbe_0 * tke_0;

            to_z_z_xz_0[k] = -2.0 * to_x_z[k] * tke_0 + 4.0 * to_xzz_z[k] * tbe_0 * tke_0;

            to_z_z_yy_0[k] = 4.0 * to_yyz_z[k] * tbe_0 * tke_0;

            to_z_z_yz_0[k] = -2.0 * to_y_z[k] * tke_0 + 4.0 * to_yzz_z[k] * tbe_0 * tke_0;

            to_z_z_zz_0[k] = -4.0 * to_z_z[k] * tke_0 + 4.0 * to_zzz_z[k] * tbe_0 * tke_0;
        }
    }
}

}  // namespace t2cgeom
