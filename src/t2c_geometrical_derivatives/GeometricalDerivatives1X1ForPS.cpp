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

#include "GeometricalDerivatives1X1ForPS.hpp"

namespace t2cgeom {  // t2cgeom namespace

auto
comp_prim_op_geom_11_ps(CSimdArray<double>&       pbuffer,
                        const size_t              idx_op_geom_101_ps,
                        const size_t              idx_op_sp,
                        const size_t              idx_op_dp,
                        const size_t              op_comps,
                        const CSimdArray<double>& factors,
                        const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : SP

        auto to_0_x = pbuffer.data(idx_op_sp + i * 3 + 0);

        auto to_0_y = pbuffer.data(idx_op_sp + i * 3 + 1);

        auto to_0_z = pbuffer.data(idx_op_sp + i * 3 + 2);

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

        // Set up 0-3 components of targeted buffer : PS

        auto to_x_x_x_0 = pbuffer.data(idx_op_geom_101_ps + 0 * op_comps * 3 + i * 3 + 0);

        auto to_x_x_y_0 = pbuffer.data(idx_op_geom_101_ps + 0 * op_comps * 3 + i * 3 + 1);

        auto to_x_x_z_0 = pbuffer.data(idx_op_geom_101_ps + 0 * op_comps * 3 + i * 3 + 2);

#pragma omp simd aligned(to_0_x, to_x_x_x_0, to_x_x_y_0, to_x_x_z_0, to_xx_x, to_xy_x, to_xz_x, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_x_0[k] = -2.0 * to_0_x[k] * tke_0 + 4.0 * to_xx_x[k] * tbe_0 * tke_0;

            to_x_x_y_0[k] = 4.0 * to_xy_x[k] * tbe_0 * tke_0;

            to_x_x_z_0[k] = 4.0 * to_xz_x[k] * tbe_0 * tke_0;
        }

        // Set up 3-6 components of targeted buffer : PS

        auto to_x_y_x_0 = pbuffer.data(idx_op_geom_101_ps + 1 * op_comps * 3 + i * 3 + 0);

        auto to_x_y_y_0 = pbuffer.data(idx_op_geom_101_ps + 1 * op_comps * 3 + i * 3 + 1);

        auto to_x_y_z_0 = pbuffer.data(idx_op_geom_101_ps + 1 * op_comps * 3 + i * 3 + 2);

#pragma omp simd aligned(to_0_y, to_x_y_x_0, to_x_y_y_0, to_x_y_z_0, to_xx_y, to_xy_y, to_xz_y, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_x_0[k] = -2.0 * to_0_y[k] * tke_0 + 4.0 * to_xx_y[k] * tbe_0 * tke_0;

            to_x_y_y_0[k] = 4.0 * to_xy_y[k] * tbe_0 * tke_0;

            to_x_y_z_0[k] = 4.0 * to_xz_y[k] * tbe_0 * tke_0;
        }

        // Set up 6-9 components of targeted buffer : PS

        auto to_x_z_x_0 = pbuffer.data(idx_op_geom_101_ps + 2 * op_comps * 3 + i * 3 + 0);

        auto to_x_z_y_0 = pbuffer.data(idx_op_geom_101_ps + 2 * op_comps * 3 + i * 3 + 1);

        auto to_x_z_z_0 = pbuffer.data(idx_op_geom_101_ps + 2 * op_comps * 3 + i * 3 + 2);

#pragma omp simd aligned(to_0_z, to_x_z_x_0, to_x_z_y_0, to_x_z_z_0, to_xx_z, to_xy_z, to_xz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_x_0[k] = -2.0 * to_0_z[k] * tke_0 + 4.0 * to_xx_z[k] * tbe_0 * tke_0;

            to_x_z_y_0[k] = 4.0 * to_xy_z[k] * tbe_0 * tke_0;

            to_x_z_z_0[k] = 4.0 * to_xz_z[k] * tbe_0 * tke_0;
        }

        // Set up 9-12 components of targeted buffer : PS

        auto to_y_x_x_0 = pbuffer.data(idx_op_geom_101_ps + 3 * op_comps * 3 + i * 3 + 0);

        auto to_y_x_y_0 = pbuffer.data(idx_op_geom_101_ps + 3 * op_comps * 3 + i * 3 + 1);

        auto to_y_x_z_0 = pbuffer.data(idx_op_geom_101_ps + 3 * op_comps * 3 + i * 3 + 2);

#pragma omp simd aligned(to_0_x, to_xy_x, to_y_x_x_0, to_y_x_y_0, to_y_x_z_0, to_yy_x, to_yz_x, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_x_0[k] = 4.0 * to_xy_x[k] * tbe_0 * tke_0;

            to_y_x_y_0[k] = -2.0 * to_0_x[k] * tke_0 + 4.0 * to_yy_x[k] * tbe_0 * tke_0;

            to_y_x_z_0[k] = 4.0 * to_yz_x[k] * tbe_0 * tke_0;
        }

        // Set up 12-15 components of targeted buffer : PS

        auto to_y_y_x_0 = pbuffer.data(idx_op_geom_101_ps + 4 * op_comps * 3 + i * 3 + 0);

        auto to_y_y_y_0 = pbuffer.data(idx_op_geom_101_ps + 4 * op_comps * 3 + i * 3 + 1);

        auto to_y_y_z_0 = pbuffer.data(idx_op_geom_101_ps + 4 * op_comps * 3 + i * 3 + 2);

#pragma omp simd aligned(to_0_y, to_xy_y, to_y_y_x_0, to_y_y_y_0, to_y_y_z_0, to_yy_y, to_yz_y, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_x_0[k] = 4.0 * to_xy_y[k] * tbe_0 * tke_0;

            to_y_y_y_0[k] = -2.0 * to_0_y[k] * tke_0 + 4.0 * to_yy_y[k] * tbe_0 * tke_0;

            to_y_y_z_0[k] = 4.0 * to_yz_y[k] * tbe_0 * tke_0;
        }

        // Set up 15-18 components of targeted buffer : PS

        auto to_y_z_x_0 = pbuffer.data(idx_op_geom_101_ps + 5 * op_comps * 3 + i * 3 + 0);

        auto to_y_z_y_0 = pbuffer.data(idx_op_geom_101_ps + 5 * op_comps * 3 + i * 3 + 1);

        auto to_y_z_z_0 = pbuffer.data(idx_op_geom_101_ps + 5 * op_comps * 3 + i * 3 + 2);

#pragma omp simd aligned(to_0_z, to_xy_z, to_y_z_x_0, to_y_z_y_0, to_y_z_z_0, to_yy_z, to_yz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_x_0[k] = 4.0 * to_xy_z[k] * tbe_0 * tke_0;

            to_y_z_y_0[k] = -2.0 * to_0_z[k] * tke_0 + 4.0 * to_yy_z[k] * tbe_0 * tke_0;

            to_y_z_z_0[k] = 4.0 * to_yz_z[k] * tbe_0 * tke_0;
        }

        // Set up 18-21 components of targeted buffer : PS

        auto to_z_x_x_0 = pbuffer.data(idx_op_geom_101_ps + 6 * op_comps * 3 + i * 3 + 0);

        auto to_z_x_y_0 = pbuffer.data(idx_op_geom_101_ps + 6 * op_comps * 3 + i * 3 + 1);

        auto to_z_x_z_0 = pbuffer.data(idx_op_geom_101_ps + 6 * op_comps * 3 + i * 3 + 2);

#pragma omp simd aligned(to_0_x, to_xz_x, to_yz_x, to_z_x_x_0, to_z_x_y_0, to_z_x_z_0, to_zz_x, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_x_0[k] = 4.0 * to_xz_x[k] * tbe_0 * tke_0;

            to_z_x_y_0[k] = 4.0 * to_yz_x[k] * tbe_0 * tke_0;

            to_z_x_z_0[k] = -2.0 * to_0_x[k] * tke_0 + 4.0 * to_zz_x[k] * tbe_0 * tke_0;
        }

        // Set up 21-24 components of targeted buffer : PS

        auto to_z_y_x_0 = pbuffer.data(idx_op_geom_101_ps + 7 * op_comps * 3 + i * 3 + 0);

        auto to_z_y_y_0 = pbuffer.data(idx_op_geom_101_ps + 7 * op_comps * 3 + i * 3 + 1);

        auto to_z_y_z_0 = pbuffer.data(idx_op_geom_101_ps + 7 * op_comps * 3 + i * 3 + 2);

#pragma omp simd aligned(to_0_y, to_xz_y, to_yz_y, to_z_y_x_0, to_z_y_y_0, to_z_y_z_0, to_zz_y, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_x_0[k] = 4.0 * to_xz_y[k] * tbe_0 * tke_0;

            to_z_y_y_0[k] = 4.0 * to_yz_y[k] * tbe_0 * tke_0;

            to_z_y_z_0[k] = -2.0 * to_0_y[k] * tke_0 + 4.0 * to_zz_y[k] * tbe_0 * tke_0;
        }

        // Set up 24-27 components of targeted buffer : PS

        auto to_z_z_x_0 = pbuffer.data(idx_op_geom_101_ps + 8 * op_comps * 3 + i * 3 + 0);

        auto to_z_z_y_0 = pbuffer.data(idx_op_geom_101_ps + 8 * op_comps * 3 + i * 3 + 1);

        auto to_z_z_z_0 = pbuffer.data(idx_op_geom_101_ps + 8 * op_comps * 3 + i * 3 + 2);

#pragma omp simd aligned(to_0_z, to_xz_z, to_yz_z, to_z_z_x_0, to_z_z_y_0, to_z_z_z_0, to_zz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_x_0[k] = 4.0 * to_xz_z[k] * tbe_0 * tke_0;

            to_z_z_y_0[k] = 4.0 * to_yz_z[k] * tbe_0 * tke_0;

            to_z_z_z_0[k] = -2.0 * to_0_z[k] * tke_0 + 4.0 * to_zz_z[k] * tbe_0 * tke_0;
        }
    }
}

}  // namespace t2cgeom
