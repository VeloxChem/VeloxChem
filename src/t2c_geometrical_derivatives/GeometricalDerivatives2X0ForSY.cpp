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

#include "GeometricalDerivatives2X0ForSY.hpp"

namespace t2cgeom {  // t2cgeom namespace

auto
comp_prim_op_geom_20_sx(CSimdArray<double>& pbuffer,
                        const size_t        idx_op_geom_200_ss,
                        const size_t        idx_op_ss,
                        const size_t        idx_op_ds,
                        const size_t        op_comps,
                        const size_t        ket_comps,
                        const double        a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    for (size_t i = 0; i < op_comps; i++)
    {
        for (size_t j = 0; j < ket_comps; j++)
        {
            // Set up components of auxiliary buffer : SS

            auto to_0_0 = pbuffer.data(idx_op_ss + i * 1 * ket_comps + 0 * ket_comps + j);

            // Set up components of auxiliary buffer : DS

            auto to_xx_0 = pbuffer.data(idx_op_ds + i * 6 * ket_comps + 0 * ket_comps + j);

            auto to_xy_0 = pbuffer.data(idx_op_ds + i * 6 * ket_comps + 1 * ket_comps + j);

            auto to_xz_0 = pbuffer.data(idx_op_ds + i * 6 * ket_comps + 2 * ket_comps + j);

            auto to_yy_0 = pbuffer.data(idx_op_ds + i * 6 * ket_comps + 3 * ket_comps + j);

            auto to_yz_0 = pbuffer.data(idx_op_ds + i * 6 * ket_comps + 4 * ket_comps + j);

            auto to_zz_0 = pbuffer.data(idx_op_ds + i * 6 * ket_comps + 5 * ket_comps + j);

            // Set up components of targeted buffer : SS

            auto to_xx_0_0_0 = pbuffer.data(idx_op_geom_200_ss + 0 * op_comps * 1 * ket_comps + i * 1 * ket_comps + 0 * ket_comps + j);

            auto to_xy_0_0_0 = pbuffer.data(idx_op_geom_200_ss + 1 * op_comps * 1 * ket_comps + i * 1 * ket_comps + 0 * ket_comps + j);

            auto to_xz_0_0_0 = pbuffer.data(idx_op_geom_200_ss + 2 * op_comps * 1 * ket_comps + i * 1 * ket_comps + 0 * ket_comps + j);

            auto to_yy_0_0_0 = pbuffer.data(idx_op_geom_200_ss + 3 * op_comps * 1 * ket_comps + i * 1 * ket_comps + 0 * ket_comps + j);

            auto to_yz_0_0_0 = pbuffer.data(idx_op_geom_200_ss + 4 * op_comps * 1 * ket_comps + i * 1 * ket_comps + 0 * ket_comps + j);

            auto to_zz_0_0_0 = pbuffer.data(idx_op_geom_200_ss + 5 * op_comps * 1 * ket_comps + i * 1 * ket_comps + 0 * ket_comps + j);

#pragma omp simd aligned(to_0_0,          \
                             to_xx_0,     \
                             to_xx_0_0_0, \
                             to_xy_0,     \
                             to_xy_0_0_0, \
                             to_xz_0,     \
                             to_xz_0_0_0, \
                             to_yy_0,     \
                             to_yy_0_0_0, \
                             to_yz_0,     \
                             to_yz_0_0_0, \
                             to_zz_0,     \
                             to_zz_0_0_0 : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                const double tbe_0 = a_exp;

                to_xx_0_0_0[k] = -2.0 * to_0_0[k] * tbe_0 + 4.0 * to_xx_0[k] * tbe_0 * tbe_0;

                to_xy_0_0_0[k] = 4.0 * to_xy_0[k] * tbe_0 * tbe_0;

                to_xz_0_0_0[k] = 4.0 * to_xz_0[k] * tbe_0 * tbe_0;

                to_yy_0_0_0[k] = -2.0 * to_0_0[k] * tbe_0 + 4.0 * to_yy_0[k] * tbe_0 * tbe_0;

                to_yz_0_0_0[k] = 4.0 * to_yz_0[k] * tbe_0 * tbe_0;

                to_zz_0_0_0[k] = -2.0 * to_0_0[k] * tbe_0 + 4.0 * to_zz_0[k] * tbe_0 * tbe_0;
            }
        }
    }
}

}  // namespace t2cgeom
