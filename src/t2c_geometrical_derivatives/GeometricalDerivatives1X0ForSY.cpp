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

#include "GeometricalDerivatives1X0ForSY.hpp"

namespace t2cgeom {  // t2cgeom namespace

auto
comp_prim_op_geom_10_sx(CSimdArray<double>& pbuffer,
                        const size_t        idx_op_geom_100_ss,
                        const size_t        idx_op_ps,
                        const size_t        op_comps,
                        const size_t        ket_comps,
                        const double        a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    for (size_t i = 0; i < op_comps; i++)
    {
        for (size_t j = 0; j < ket_comps; j++)
        {
            // Set up components of auxiliary buffer : PS

            auto to_x_0 = pbuffer.data(idx_op_ps + i * 3 * ket_comps + 0 * ket_comps + j);

            auto to_y_0 = pbuffer.data(idx_op_ps + i * 3 * ket_comps + 1 * ket_comps + j);

            auto to_z_0 = pbuffer.data(idx_op_ps + i * 3 * ket_comps + 2 * ket_comps + j);

            // Set up components of targeted buffer : SS

            auto to_x_0_0_0 = pbuffer.data(idx_op_geom_100_ss + 0 * op_comps * 1 * ket_comps + i * 1 * ket_comps + 0 * ket_comps + j);

            auto to_y_0_0_0 = pbuffer.data(idx_op_geom_100_ss + 1 * op_comps * 1 * ket_comps + i * 1 * ket_comps + 0 * ket_comps + j);

            auto to_z_0_0_0 = pbuffer.data(idx_op_geom_100_ss + 2 * op_comps * 1 * ket_comps + i * 1 * ket_comps + 0 * ket_comps + j);

#pragma omp simd aligned(to_x_0, to_x_0_0_0, to_y_0, to_y_0_0_0, to_z_0, to_z_0_0_0 : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                const double tbe_0 = a_exp;

                to_x_0_0_0[k] = 2.0 * to_x_0[k] * tbe_0;

                to_y_0_0_0[k] = 2.0 * to_y_0[k] * tbe_0;

                to_z_0_0_0[k] = 2.0 * to_z_0[k] * tbe_0;
            }
        }
    }
}

}  // namespace t2cgeom
