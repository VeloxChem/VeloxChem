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

#include "ElectricDipoleMomentumPrimRecSS.hpp"

namespace diprec {  // diprec namespace

auto
comp_prim_electric_dipole_momentum_ss(CSimdArray<double>&       pbuffer,
                                      const size_t              idx_dip_ss,
                                      const size_t              idx_ovl_ss,
                                      const CSimdArray<double>& factors,
                                      const size_t              idx_rpc) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(PC) distances

    auto pc_x = factors.data(idx_rpc);

    auto pc_y = factors.data(idx_rpc + 1);

    auto pc_z = factors.data(idx_rpc + 2);

    /// Set up components of auxiliary buffer : SS

    auto ts_0_0 = pbuffer.data(idx_ovl_ss);

    /// Set up components of targeted buffer : SS

    auto tr_x_0_0 = pbuffer.data(idx_dip_ss);

    auto tr_y_0_0 = pbuffer.data(idx_dip_ss + 1);

    auto tr_z_0_0 = pbuffer.data(idx_dip_ss + 2);

#pragma omp simd aligned(tr_x_0_0, tr_y_0_0, tr_z_0_0, ts_0_0, pc_x, pc_y, pc_z : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fact = ts_0_0[i];

        tr_x_0_0[i] = fact * pc_x[i];

        tr_y_0_0[i] = fact * pc_y[i];

        tr_z_0_0[i] = fact * pc_z[i];
    }
}

}  // namespace diprec
