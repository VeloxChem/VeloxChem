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

#include "GeometricalDerivatives1X1ForSS.hpp"

namespace t2cgeom {  // t2cgeom namespace

auto
comp_prim_op_geom_11_ss(CSimdArray<double>&       pbuffer,
                        const size_t              idx_op_geom_101_ss,
                        const size_t              idx_op_pp,
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

        // Set up components of targeted buffer : SS

        auto to_x_x_0_0 = pbuffer.data(idx_op_geom_101_ss + 0 * op_comps * 1 + i * 1 + 0);

        auto to_x_y_0_0 = pbuffer.data(idx_op_geom_101_ss + 1 * op_comps * 1 + i * 1 + 0);

        auto to_x_z_0_0 = pbuffer.data(idx_op_geom_101_ss + 2 * op_comps * 1 + i * 1 + 0);

        auto to_y_x_0_0 = pbuffer.data(idx_op_geom_101_ss + 3 * op_comps * 1 + i * 1 + 0);

        auto to_y_y_0_0 = pbuffer.data(idx_op_geom_101_ss + 4 * op_comps * 1 + i * 1 + 0);

        auto to_y_z_0_0 = pbuffer.data(idx_op_geom_101_ss + 5 * op_comps * 1 + i * 1 + 0);

        auto to_z_x_0_0 = pbuffer.data(idx_op_geom_101_ss + 6 * op_comps * 1 + i * 1 + 0);

        auto to_z_y_0_0 = pbuffer.data(idx_op_geom_101_ss + 7 * op_comps * 1 + i * 1 + 0);

        auto to_z_z_0_0 = pbuffer.data(idx_op_geom_101_ss + 8 * op_comps * 1 + i * 1 + 0);

#pragma omp simd aligned(to_x_x,         \
                             to_x_x_0_0, \
                             to_x_y,     \
                             to_x_y_0_0, \
                             to_x_z,     \
                             to_x_z_0_0, \
                             to_y_x,     \
                             to_y_x_0_0, \
                             to_y_y,     \
                             to_y_y_0_0, \
                             to_y_z,     \
                             to_y_z_0_0, \
                             to_z_x,     \
                             to_z_x_0_0, \
                             to_z_y,     \
                             to_z_y_0_0, \
                             to_z_z,     \
                             to_z_z_0_0, \
                             b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_0_0[k] = 4.0 * to_x_x[k] * tbe_0 * tke_0;

            to_x_y_0_0[k] = 4.0 * to_x_y[k] * tbe_0 * tke_0;

            to_x_z_0_0[k] = 4.0 * to_x_z[k] * tbe_0 * tke_0;

            to_y_x_0_0[k] = 4.0 * to_y_x[k] * tbe_0 * tke_0;

            to_y_y_0_0[k] = 4.0 * to_y_y[k] * tbe_0 * tke_0;

            to_y_z_0_0[k] = 4.0 * to_y_z[k] * tbe_0 * tke_0;

            to_z_x_0_0[k] = 4.0 * to_z_x[k] * tbe_0 * tke_0;

            to_z_y_0_0[k] = 4.0 * to_z_y[k] * tbe_0 * tke_0;

            to_z_z_0_0[k] = 4.0 * to_z_z[k] * tbe_0 * tke_0;
        }
    }
}

}  // namespace t2cgeom
