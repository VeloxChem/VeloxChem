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

#include "ThreeCenterElectronRepulsionPrimRecSSS.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_sss(CSimdArray<double>&       pbuffer,
                                 const size_t              idx_eri_0_sss,
                                 CSimdArray<double>&       factors,
                                 const size_t              idx_ovl,
                                 const CSimdArray<double>& bf_data,
                                 const size_t              idx_bvals) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    /// Set up components of targeted buffer : prim_buffer_0_sss

    auto g_0_0_0_0 = pbuffer.data(idx_eri_0_sss);

    /// Set up overlap factors

    auto fovl_acd = factors.data(idx_ovl);

    /// Set up Boys function values

    auto bf_values = bf_data.data(idx_bvals);

#pragma omp simd aligned(g_0_0_0_0, fovl_acd, bf_values : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_0_0_0_0[i] = fovl_acd[i] * bf_values[i];
    }
    
}

} // t3ceri namespace

