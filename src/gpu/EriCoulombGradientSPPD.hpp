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

#ifndef EriCoulombGradientSPPD_hpp
#define EriCoulombGradientSPPD_hpp

#include <cstdint>

#include "GpuConstants.hpp"

namespace gpu {  // gpu namespace


__global__ void __launch_bounds__(TILE_SIZE_J)
computeCoulombGradientSPPD_I_0(double*         grad_x,
                       const uint32_t  grad_cart_ind,
                       const double    prefac_coulomb,
                       const double*   s_prim_info,
                       const uint32_t  s_prim_count,
                       const double*   p_prim_info,
                       const uint32_t  p_prim_count,
                       const double*   d_prim_info,
                       const uint32_t  d_prim_count,
                       const double*   sp_mat_D_local,
                       const double*   pd_mat_D,
                       const double*   sp_mat_Q_local,
                       const double*   pd_mat_Q,
                       const uint32_t* sp_first_inds_local,
                       const uint32_t* sp_second_inds_local,
                       const double*   sp_pair_data_local,
                       const uint32_t  sp_prim_pair_count_local,
                       const uint32_t* pd_first_inds,
                       const uint32_t* pd_second_inds,
                       const double*   pd_pair_data,
                       const uint32_t  pd_prim_pair_count,
                       const uint32_t* prim_cart_ao_to_atom_inds,
                       const double*   boys_func_table,
                       const double*   boys_func_ft,
                       const double    eri_threshold);

__global__ void __launch_bounds__(TILE_SIZE_J)
computeCoulombGradientSPPD_J_0(double*         grad_x,
                       const uint32_t  grad_cart_ind,
                       const double    prefac_coulomb,
                       const double*   s_prim_info,
                       const uint32_t  s_prim_count,
                       const double*   p_prim_info,
                       const uint32_t  p_prim_count,
                       const double*   d_prim_info,
                       const uint32_t  d_prim_count,
                       const double*   sp_mat_D_local,
                       const double*   pd_mat_D,
                       const double*   sp_mat_Q_local,
                       const double*   pd_mat_Q,
                       const uint32_t* sp_first_inds_local,
                       const uint32_t* sp_second_inds_local,
                       const double*   sp_pair_data_local,
                       const uint32_t  sp_prim_pair_count_local,
                       const uint32_t* pd_first_inds,
                       const uint32_t* pd_second_inds,
                       const double*   pd_pair_data,
                       const uint32_t  pd_prim_pair_count,
                       const uint32_t* prim_cart_ao_to_atom_inds,
                       const double*   boys_func_table,
                       const double*   boys_func_ft,
                       const double    eri_threshold);

}  // namespace gpu

#endif
